import Analysis
import cPickle as pickle
import Tools
import multiprocessing as mp 
import pyfits
import numpy as np
import h5py
import sys


class GlobalFitter():

    def __init__(self):
        pass

    def GenDiffuse(self, analysis, basedir, tag, verbosity=0, multiplier=1., nrings=9,
                    E_subsample=1, fixSpectrum=True, fix_xco=False):
            """
            This method adds all the galprop templates to the analysis. 
      
            :param analysis: The analysis object being used. 
            :param basedir: Directory containing the galprop output. 
            :param tag: Tag for the galprop file.  This is the part between '_54_' and '.gz' of the galprop output. 
            :param verbosity: 0 is quiet, >1 prints status.
            :param multiplier: PSF is a gaussian with sigma=r_68*multiplier, where r_68 is the 68% containment radius.  
            :param nrings: Number of rings in galprop gas map.  Normally 9 for reduced ring models included in the support_material. 
            :param fix_xco = False:  If True, the XCO rings are not fit, but instead are fixed to their galprop values.
            :param E_subsample=1: can subsample the flux density over bin.  Not really necessary to be higher than 1. 
            :param fixSpectrum=True: Fix the spectrum of all components to the galprop output. 
            """

            if verbosity>0:
                print 'Loading FITS'
            # Load the energies from the galprop output. 
            energies = pyfits.open(basedir+'/bremss_healpix_54_'+tag+'.gz')[2].data.field(0)
            comps, comps_new = {}, {}
    #         # For some reason, older versions of galprop files have slightly different data structures.  This try/except
    #         # will detect the right one to use. 
             
    #         try:
    #             comps['ics'] = pyfits.open(basedir+'/ics_isotropic_healpix_54_'+tag+'.gz')[1].data.field(0).T
    #             nside_in = np.sqrt(comps['ics'].shape[1]/12)
    #             comps['pi0'] = pyfits.open(basedir+'/pi0_decay_healpix_54_'+tag+'.gz')[1].data.field(0).T
    #             comps['brem'] = pyfits.open(basedir+'/bremss_healpix_54_'+tag+'.gz')[1].data.field(0).T

            #except:
            # Convert the galprop fits files into a usable numpy array. 
            def ReadFits(fname, length):
                d = pyfits.open(fname)[1].data
                return np.array([d.field(i) for i in range(length)])
            
            # Add up the HI and HII rings into a single template since nothing there is varying.
            pi0HIHII = np.zeros((len(energies), 12*analysis.nside**2))
            bremHIHII = np.zeros((len(energies), 12*analysis.nside**2))

            for i_ring in range(1,nrings+1):
                print "Adding HI/HII ring", i_ring
                bremHIHII += ReadFits(basedir+'/bremss_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                bremHIHII += ReadFits(basedir+'/bremss_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                pi0HIHII += ReadFits(basedir+'/pi0_decay_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                pi0HIHII += ReadFits(basedir+'/pi0_decay_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                if fix_xco:
                    bremHIHII += ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                    pi0HIHII += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))

            # By default we merge bremsstrahlung with the pi0 at a ratio 1.25 to 1. 
            # Otherwise the XCO fits are highly degenerate. 
            comps['pi0HIHII'] = pi0HIHII + 1.25*bremHIHII
            # comps['pi0HIHII'] = pi0HIHII
            # comps['bremHIHII'] = bremHIHII
            comps_new['pi0HIHII'] =  np.zeros((analysis.n_bins, 12*analysis.nside**2))
            # comps_new['bremHIHII'] =  np.zeros((analysis.n_bins, 12*analysis.nside**2))
            
            # Loop over rings and add the H2 contributions from pi0 and brems.
            for i_ring in range(1,nrings+1):
                if not fix_xco:
                    print "Adding H2 ring", i_ring
                    brem = ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                    pi = ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies))
                    comps['pi0_H2_'+str(i_ring)] = pi + 1.25*brem

                    # comps['pi0_H2_'+str(i_ring)] = pi 
                    # comps['brem_H2_'+str(i_ring)]= brem
                    comps_new['pi0_H2_'+str(i_ring)] =  np.zeros((analysis.n_bins, 12*analysis.nside**2))
                    # comps_new['brem_H2_'+str(i_ring)] = np.zeros((analysis.n_bins, 12*analysis.nside**2))
            
            # Add the ICS components.
            comps['ics'] = np.zeros((len(energies), 12*analysis.nside**2))
            comps_new['ics'] = np.zeros((analysis.n_bins, 12*analysis.nside**2))
            for i_ics in range(1,4):
                print "Adding ics", i_ics
                comps['ics'] += ReadFits(basedir+'/ics_isotropic_comp_'+str(i_ics)+'_healpix_54_'+tag+'.gz', len(energies))
    #             comps['ics_'+str(i_ics)] = ReadFits(basedir+'/ics_isotropic_comp_'+str(i_ics)+'_healpix_54_'+tag+'.gz', len(energies))
    #             comps_new['ics_'+str(i_ics)] = np.zeros((analysis.n_bins, 12*analysis.nside**2))
            
            # Determine resolution of input maps. 
            nside_in = np.sqrt(comps['pi0HIHII'].shape[1]/12)
            
            #---------------------------------------------------------------------------------
            # Now we integrate each model over the energy bins...
            #
            # Multiprocessing for speed. There is an async callback which applies each result to
            # the arrays.  Not sure why RunAsync needs new thread pool for each component, but this
            # works and decreases memory footprint.
            def callback(result):
                idx, comp, dat = result
                comps_new[comp][idx] = dat

            def RunAsync(component):
                p = mp.Pool(mp.cpu_count())
                for i_E in range(analysis.n_bins):
                    p.apply_async(Tools.AsyncInterpolateHealpix,
                                  [comps[component], energies, analysis.bin_edges[i_E], analysis.bin_edges[i_E+1],
                                   i_E, component, E_subsample, analysis.nside],
                                  callback=callback)
                p.close()
                p.join()

            # For each component, run the async sampling/sizing.  
            # This converts the galprop maps into the correct spectral binning and spatial resolution.
            # 
            for key in comps:
                if verbosity>0:
                    print 'Integrating and Resampling', key, 'templates...'
                    sys.stdout.flush()
                for i_E in range(analysis.n_bins):
                    comps_new[key][i_E] = Tools.InterpolateHealpix(comps[key], energies,  
                        analysis.bin_edges[i_E], analysis.bin_edges[i_E+1], E_bins=E_subsample, nside_out=analysis.nside)
                # Parallel version. (very memory hungry)
                #RunAsync(key)


            #---------------------------------------------------------------------------------
            # Now we just need to add the templates to the active template stack
            print 'Adding Templates to stack'
            
            analysis.AddTemplate(name='pi0HIHII', healpixCube=comps_new['pi0HIHII'].clip(0), fixSpectrum=fixSpectrum, fixNorm=False,
                               value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)
            # analysis.AddTemplate(name='bremHIHII', healpixCube=comps_new['bremHIHII'], fixSpectrum=fixSpectrum, fixNorm=False,
                               # value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)
            
            for i_ring in range(1,nrings+1):
                if not fix_xco:
                    analysis.AddTemplate(name='pi0_H2_'+str(i_ring), healpixCube=comps_new['pi0_H2_'+str(i_ring)].clip(0), fixSpectrum=fixSpectrum, fixNorm=False,
                                   value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)
                    
                    # analysis.AddTemplate(name='brem_H2_'+str(i_ring), healpixCube=comps_new['brem_H2_'+str(i_ring)], fixSpectrum=fixSpectrum, fixNorm=False,
                    #                value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)
                
        #         for i_ics in range(1,4):
        #             analysis.AddTemplate(name='ics_'+str(i_ics), healpixCube=comps_new['ics_'+str(i_ics)], fixSpectrum=fixSpectrum, fixNorm=False,
        #                            value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)
            analysis.AddTemplate(name='ics', healpixCube=comps_new['ics'].clip(0), fixSpectrum=fixSpectrum, fixNorm=False,
                           value=1., ApplyIRF=True,noPSF=True, sourceClass='GEN', limits=[None, None], multiplier=multiplier)

            return analysis



    def RunFit(self, A, nrings=9, limit_inner=None, fix_xco=False, tag='', output_loglike=False):
        '''Run the three stage global fit described in arXiv: 1603.06584 to fit the XCO profile in the Galaxy. 
        This expects a 9 ring galprop gas model, and must be modified if fitting another model. 
        
        :param A: Analysis object containing the model to be fit. 
        :param nrings: Number of galprop gas rings.  Currently, the fitter is limited to 9 ring model without modification.
        :param limit_inner: Place a lower limit on XCO across all rings.  Do not use without looking at code below. 
        :param fix_xco: Lock XCO to the galprop model. The fit is still performed with other components. 
        :param tag: Tag for the galprop file.  This is the part between '_54_' and '.gz' of the galprop output. 
        :param output_loglike: if True, outputs a pixel-by-pixel log-likelihood.  
        '''

        #-----------------------------------------------------
        # Templates are now added so we fit X_CO
        
        import GammaLikelihood as like
        
        if output_loglike:
            print 'loglike output ', 'loglike_'+tag+'.npy' 

        fval, res = [], []
        for key, t in A.templateList.items():
            # local 
            if key not in ['pi0_H2_7', 'Isotropic','pi0HIHII', 'ics', 'Bubbles']:
                t.fixNorm = True
                t.fixSpectrum= True
                t.limits = [0.0,200.]
                t.value=1.
            else: 
                t.fixNorm = False
                t.limits=[0., None]

                
        print 'Running Local Ring Fit...'
        A.GenSquareMask(l_range=[-180.,180], b_range=[-85.,85.], plane_mask=8.)
        m, R = A.RunLikelihood( print_level=1, precision=None, tol=1e3)[:2]
        fval.append(m.fval)
        res.append(R)


        # output loglikeliehood.
        if output_loglike:
            model = np.zeros(np.shape(A.binned_data))
            for name, t in A.templateList.items():
                model += t.healpixCube*t.value
            loglike = A.psc_weights*A.mask*(model-A.binned_data*np.log(model))


        print 'isotropic value:', A.templateList['Isotropic'].value
        
        if fix_xco is False:
            vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
            print "X_CO adjustment (this is not XCO value, it is multiplier for galdef values of xco):", vals

        #-------------------------------------------------------------------
        # Now we have fit the local X_CO (fixed).  Next we fit the outer galaxy
        if fix_xco is False:
            A.templateList['pi0_H2_7'].fixNorm = True
        A.templateList['Isotropic'].fixNorm = True
        A.templateList['Bubbles'].fixNorm = True
        # Let the outer two rings float
        if fix_xco is False:
            A.templateList['pi0_H2_8'].fixNorm = False
            A.templateList['pi0_H2_9'].fixNorm = False

        print 'Running Outer Rings Fit...'
        A.GenSquareMask(l_range=[-180,-80], b_range=[-8.,8.], plane_mask=0)
        A.GenSquareMask(l_range=[80,180], b_range=[-8.,8.], plane_mask=0, merge=True)

        m, R = A.RunLikelihood( print_level=1, precision=None, tol=1e3)[:2]
        fval.append(m.fval)
        res.append(R)

        if fix_xco is False:
            vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
            print "X_CO adjustment (this is not XCO value, it is multiplier for galdef values of xco):", vals

        # output loglikeliehood.
        if output_loglike:
            model = np.zeros(np.shape(A.binned_data))
            for name, t in A.templateList.items():
                model += t.healpixCube*t.value
            loglike += A.psc_weights*A.mask*(model-A.binned_data*np.log(model))


        #-------------------------------------------------------------------
        # Now we fit the inner galaxy X_CO.
        if fix_xco is False:
            A.templateList['pi0_H2_8'].fixNorm = True
            A.templateList['pi0_H2_9'].fixNorm = True
        
        if limit_inner is not None:
            if fix_xco is False:
                for j in range(1,10):
                    X_CO_values_GAL = [1e+19, 8.42e+19, 1.61e+20, 1.73e+20, 1.72e+20, 1.74e+20, 8.61e+19, 4.29e+20, 2.01e+21]
                    A.templateList['pi0_H2_%i'%j].limits=[limit_inner*1e19/X_CO_values_GAL[j-1],None] # set inner ring limit to not fall below inner_limit
                    #A.templateList['pi0_H2_%i'%j].limits=[limit_inner*1e19/4e19,None] # set inner ring limit to not fall below inner_limit

        # Let the inner 6 rings float
        for i in range(1,7):
            if fix_xco is False:
                A.templateList['pi0_H2_' + str(i)].fixNorm=False

        print 'Running Inner Rings Fit...'
        A.GenSquareMask(l_range=[-80,80], b_range=[-8.,8.], plane_mask=0)
        m, R = A.RunLikelihood( print_level=1, precision=None, tol=5e3)[:2]
        fval.append(m.fval)
        res.append(R)

        if fix_xco is False:
            vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
            print "X_CO adjustment (this is not XCO value, it is multiplier for galdef values of xco):", vals

        if output_loglike:
            model = np.zeros(np.shape(A.binned_data))
            for name, t in A.templateList.items():
                model += t.healpixCube*t.value
            loglike += A.psc_weights*A.mask*(model-A.binned_data*np.log(model))
            
            np.save('/pfs/carlson/GCE_sys/loglike_'+tag+'.npy', loglike)

        return m, fval, res





    def WriteHDF5(self, fname, basedir, tag, m , nrings=9, fix_xco=False):
        """
        Build the diffuse model according to the best fit parameters and writes model+metadata to an HDF5 file.
        
        :param fname: Output filename for the HDF5 file
        :param basedir: Directory with galprop output files
        :param tag: galprop tag
        :param m: iminuit object (usually just analysis.m). 
        :param nrings: number of galprop rings.
        :param fix_xco: Lock XCO to the galprop model. The fit is still performed with other components. 
        """
        modf = h5py.File(fname, 'w')

        if m is not None:
            if fix_xco is False:
                X_CO = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,nrings+1)])
            else: 
                X_CO = np.ones(nrings)
        #modf = h5py.File(fname, 'w')
        #try:
        # Generate Groups
        # template_group = modf.create_group("templates")
        # fit_group = modf.create_group("fit_results")


        # Get data dimensions
        tmp = pyfits.open(basedir+'/bremss_healpix_54_'+tag+'.gz')
        energies = tmp[2].data.field(0)
        tShape = (len(energies), tmp[1].data.shape[0])
        print tShape
        del tmp # free memory

        pi0     = modf.create_dataset("/templates/pi0", tShape, dtype='float32',compression="gzip")
        pi0_0   = modf.create_dataset("/templates/pi0_0", tShape, dtype='float32',compression="gzip")
        pi0_1   = modf.create_dataset("/templates/pi0_1", tShape, dtype='float32',compression="gzip")
        brem   = modf.create_dataset("/templates/brem", tShape, dtype='float32',compression="gzip")
        brem_0 = modf.create_dataset("/templates/brem_0", tShape, dtype='float32',compression="gzip")
        brem_1 = modf.create_dataset("/templates/brem_1", tShape, dtype='float32',compression="gzip")
        ics_opt = modf.create_dataset("/templates/ics_opt", tShape, dtype='float32',compression="gzip")
        ics_fir = modf.create_dataset("/templates/ics_fir", tShape, dtype='float32',compression="gzip")
        ics_cmb = modf.create_dataset("/templates/ics_cmb", tShape, dtype='float32',compression="gzip")
        modf.create_dataset("/templates/energies", data=energies, dtype='float32',compression="gzip")
        # Now fill in the templates one by one.
        # Add fit metadata.
        # Add galdef metadata.


        #---------------------------------------------------------------
        # Create Diffuse Template from fitting results.
        def ReadFits(fname, length):
            d = pyfits.open(fname)[1].data
            return np.array([d.field(i) for i in range(length)])

        if m is not None:
            for i_ring in range(1,nrings+1):
                print "Adding HI/HII ring", i_ring

                pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                if fix_xco is False:
                    pi0[...] += X_CO[i_ring-1]*ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                else:
                    pi0[...] += m.values['pi0HIHII']*ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

                brem[...] += m.values['pi0HIHII']*1.25*ReadFits(basedir+'/bremss_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                brem[...] += m.values['pi0HIHII']*1.25*ReadFits(basedir+'/bremss_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                if fix_xco is False:
                    brem[...] += 1.25*X_CO[i_ring-1]*ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                else: 
                    brem[...] += m.values['pi0HIHII']*1.25*X_CO[i_ring-1]*ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                
                if i_ring == 1 :
                    pi0_0[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                    brem_0[...] += ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                if i_ring == 2:
                    pi0_1[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                    brem_1[...] +=  ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

            print 'Adding ICS'
            ics_opt[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_1_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            ics_fir[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_2_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            ics_cmb[...] += m.values['ics']*ReadFits(basedir+'/ics_isotropic_comp_3_healpix_54_'+tag+'.gz', len(energies)).clip(0)

        else:
            for i_ring in range(1,nrings+1):
                print "Adding HI/HII/H2 ring", i_ring

                pi0[...] += ReadFits(basedir+'/pi0_decay_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                pi0[...] += ReadFits(basedir+'/pi0_decay_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                pi0[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

                brem[...] += ReadFits(basedir+'/bremss_HIR_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                brem[...] += ReadFits(basedir+'/bremss_HII_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                brem[...] += ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

                if i_ring == 1 :
                    pi0_0[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                    brem_0[...] += ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                if i_ring == 2:
                    pi0_1[...] += ReadFits(basedir+'/pi0_decay_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)
                    brem_1[...] +=  ReadFits(basedir+'/bremss_H2R_ring_'+str(i_ring)+'_healpix_54_'+tag+'.gz', len(energies)).clip(0)

            print 'Adding ICS'
            ics_opt[...] += ReadFits(basedir+'/ics_isotropic_comp_1_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            ics_fir[...] += ReadFits(basedir+'/ics_isotropic_comp_2_healpix_54_'+tag+'.gz', len(energies)).clip(0)
            ics_cmb[...] += ReadFits(basedir+'/ics_isotropic_comp_3_healpix_54_'+tag+'.gz', len(energies)).clip(0)

    #     except:
    #         modf.close()
        
        try: 
            modf.close()
            print 'Closed HDF5 file.' 
        except: 
            print 'Failed to close HDF5 file.' 
            pass
        return


    def AddMetadata(self, fname, basedir, tag, A, m, fval=None, fix_xco=False):
        """
        Parses the metadata from the galdef file and adds it to the HDF5 file.
        
        :param fname: Output filename for the HDF5 file
        :param basedir: Directory with Galdef files
        :param tag: galprop tag
        :param A: Analysis object
        :param m: iMinuit object containing fit.  Usually just A.m

        :param fix_xco: Lock XCO to the galprop model. The fit is still performed with other components. 
        """
        # Parse the galprop file into a dict.
        
        print 'Adding fit metadata'
        galdef_dict = {}
        with open(basedir + '/galdef_54_'+tag) as galdef:
            for line in galdef: 
                if line[0] != "#" and line.strip()!='':
                    s = line.strip('\n').split('=')
                    if len(s)<2: 
                        continue
                    key = s[0].strip()
                    if key in ["Title", 'X_CO_values', 'n_X_CO_values', 'X_CO_radius']: continue
                    galdef_dict[key] = s[1].strip().split(" ")[0]
        try: h5.close() 
        except: pass
        h5 = h5py.File(fname)
        try:
            galdef_group = h5.create_group("/galdef")
        except: 
            galdef_group = h5['/galdef']
        for key, val in galdef_dict.items():
            #print key, val
            galdef_group.attrs.create(key,val)
        
        if m is not None and A is not None:
            try:
                fit_results = h5.create_group("/fit_results/global")
            except: 
                fit_results = h5['/fit_results/global']
            

            

            if fix_xco is False:
                vals = np.array([m.values['pi0_H2_'+str(i)] for i in range(1,10)])
                #vals_HIH2 = np.array([m.values['pi0HIHII_'+str(i)] for i in range(1,10)])
            else: 
                vals = np.ones(9)

            fa = fit_results.attrs
            fa.create('globalvalues', m.values.items())
            fa.create('globalvaluesUnc', m.errors.items())
            fa.create('globalfval', m.fval)

            if fval is not None:
                fa.create('localfval', fval[0])            
                fa.create('outerfval', fval[1])            
                fa.create('innerfval', fval[2])            

            fa.create('global_XCO', vals)
            #fa.create('global_XHIHII', vals_HIH2)
            fa.create('globale_bins', A.bin_edges)
            fa.create('globalirf', A.irf)
            fa.create('globalevclass', A.evclass)
            if A.convtype is None:
                A.convtype = -1
            fa.create('globalconvtype', A.convtype)
            fa.create('globalphfile', A.phfile)
            fa.create('globaltag', A.tag)
            h5.create_dataset('/fit_results/globalmask', data=A.mask, dtype='float32')
        h5.close()
    