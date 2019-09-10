{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Level-1 Pipeline using MIRI data from MAST"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Table of Contents:\n",
    "> * [Detector1 Pipeline](#detector1_pipeline)\n",
    "> * [Resources and Documentation](#resources)\n",
    "> * [Download Data from MAST](#download_from_mast)\n",
    "> * [Run Pipeline with Default Configuration](#pipeline_with_defaults)\n",
    "> * [About Configuration Files](#pipeline_configs)\n",
    "> * [Run Pipeline with Configuration Files](#pipeline_with_cfgs)\n",
    "> * [Run Pipeline with Parameters Set Programmatically](#pipeline_no_configs)\n",
    "> * [Run Individual Steps with Configuration Files](#steps_with_config_files)\n",
    "> * [Run Individual Steps with Parameters Set Programmatically](#steps_no_configs)\n",
    "\n",
    "***\n",
    "<a id=detector1_pipeline></a>\n",
    "##  Detector1 Pipeline\n",
    "\n",
    "Stage 1 consists of detector-level corrections that are performed on a group-by-group basis, followed by ramp fitting. \n",
    "\n",
    "More information can be found at: https://jwst-pipeline.readthedocs.io/en/latest/jwst/pipeline/calwebb_detector1.html#calwebb-detector1\n",
    "\n",
    "> **Inputs**: The inputs to stage 1 processing will usually be level-1b raw files.\n",
    "\n",
    "> **Outputs**: The output of stage 1 processing is a countrate image per exposure, or per integration for some modes.\n",
    "    \n",
    "### Level 1 pipeline:\n",
    "\n",
    "**Calwebb Detector1**   (jwst.pipeline, calwebb_detector1, Detector1Pipeline)   (calwebb_detector1.cfg)\n",
    "\n",
    "### Level 1 pipeline steps:\n",
    "\n",
    "**Group Scale** (jwst.group_scale, group_scale_step, GroupScaleStep)   (group_scale.cfg)\n",
    "\n",
    "**DQInit** (jwst.dq_init, dq_init_step, DQInitStep)    (dq_init.cfg)\n",
    "\n",
    "**Saturation** (jwst.saturation, saturation_step, SaturationStep)    (saturation.cfg)\n",
    "\n",
    "**IPC**   (jwst.ipc,  ipc_step,  IPCStep)    (ipc.cfg)\n",
    "\n",
    "**Super Bias**   (jwst.superbias,  superbias_step,  SuperBiasStep)    (superbias.cfg)\n",
    "\n",
    "**Refpix** (jwst.refpix,  refpix_step,  RefpixStep)    (refpix.cfg)\n",
    "\n",
    "**RSCD**   (jwst.rscd,  rscd_step,  RSCD_Step)   (rscd.cfd)\n",
    "\n",
    "**First Frame**  (jwst.firstframe, firstframe_step,  FirstFrameStep)   (firstframe.cfg)\n",
    "\n",
    "**Last Frame**  (jwst.lastframe,  lastframe_step, LastFrameStep)   (lastframe.cfg)\n",
    "\n",
    "**Linearity**  (jwst.linearity,  linearity_step, LinearityStep)   (linearity.cfg)\n",
    "\n",
    "**Dark Current**  (jwst.dark_current,  dark_current_step,  DarkCurrentStep)   (dark_current.cfg)\n",
    "\n",
    "**Peristence**  (jwst.persistence,  persistence_step, PersistenceStep)   (persistence.cfg)\n",
    "\n",
    "**Jump**  (jwst.jump, jump_step,  JumpStep)  (jump_step.cfg)\n",
    "\n",
    "**RampFit**   (jwst.ramp_fit,  ramp_fit_step,  RampFitStep)   (ramp_fit_step.cfg)\n",
    "\n",
    "**Gain Scale**  (jwst.gain_scale,  gain_scale_step,  GainScaleStep)   (gain_scale_step.cfg)\n",
    "\n",
    "(for more information on individual steps see: https://jwst-pipeline.readthedocs.io/en/latest/jwst/package_index.html)\n",
    "\n",
    "***\n",
    "<a id='resources'></a>\n",
    "## Resources and Documentation\n",
    "\n",
    "There are several different places to find information on installing and running the pipeline. This notebook will give a shortened description of the steps pulled from the detailed pipeline information pages, but to find more in-depth instructions use the links below. \n",
    "\n",
    ">1. JDox: https://jwst-docs.stsci.edu/display/JDAT/JWST+Data+Reduction+Pipeline\n",
    ">2. Installation page: http://astroconda.readthedocs.io/en/latest/releases.html#pipeline-install\n",
    ">3. Detailed pipeline information: https://jwst-pipeline.readthedocs.io/en/latest/jwst/introduction.html\n",
    ">4. Help Desk (click on Pipeline Support): https://stsci.service-now.com/jwst?id=sc_category\n",
    ">5. GitHub README installation instructions: https://github.com/spacetelescope/jwst/blob/master/README.md\n",
    "\n",
    "\n",
    "If this is your first time trying to run the pipeline from a jupyter notebook, you need to install the jupyter notebook in your pipeline environment:\n",
    ">1. In a new terminal, change the directory to your working directory, terminal command: cd [your working directory]\n",
    ">2. Terminal command: source activate jwst_dev\n",
    "(or whatever your environment name for the pipeline is)\n",
    ">3. Terminal command: conda install jupyter\n",
    ">4. Terminal command: jupyter notebook\n",
    "\n",
    "**NOTE:** During your first run CRDS may download and cache reference files in $HOME/crds_cache.  On subsequent runs cached files will be used."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Imports"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "WARNING: AstropyDeprecationWarning: astropy.extern.six will be removed in 4.0, use the six module directly if it is still needed [astropy.extern.six]\n"
     ]
    }
   ],
   "source": [
    "from astroquery.mast import Observations\n",
    "from jwst.pipeline import Detector1Pipeline\n",
    "from astropy.io import fits\n",
    "import matplotlib.pyplot as plt\n",
    "%matplotlib inline"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Quick & Dirty Plots"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 52,
   "metadata": {},
   "outputs": [],
   "source": [
    "def plot_image(obj, index=None):\n",
    "    if isinstance(obj, str):\n",
    "        data = fits.getdata(obj)\n",
    "    else:\n",
    "        data = obj.data\n",
    "    image = data[index] if index is not None else data\n",
    "    fig = plt.figure(figsize=(12,12))\n",
    "    ax = fig.add_subplot(1, 1, 1)\n",
    "    im = ax.imshow(image, origin='lower') # , norm=norm)\n",
    "    fig.colorbar(im)\n",
    "    plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id=\"download_from_mast\"></a>\n",
    "## Download Data from MAST\n",
    "\n",
    "Test data can be downloaded from MAST using Astroquery"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 35,
   "metadata": {},
   "outputs": [],
   "source": [
    "DATA_URL = \"mast:JWST/product/jw00817007001_02104_00001_nrs1_uncal.fits\""
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 36,
   "metadata": {},
   "outputs": [],
   "source": [
    "def download_mast(data_url, target):\n",
    "    Observations._download_file(\"https://pwjwdmsauiweb.stsci.edu/portal/Download/file\" + \"?uri=\" + data_url, target)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 37,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Downloading URL https://pwjwdmsauiweb.stsci.edu/portal/Download/file?uri=mast:JWST/product/jw00817007001_02104_00001_nrs1_uncal.fits to test.fits ... [Done]\n"
     ]
    }
   ],
   "source": [
    "download_mast(DATA_URL, \"test.fits\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The raw images are 2k x 2k detector images in 16-bit counts."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 48,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 864x864 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_image(\"test.fits\", index=(0,0))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id=\"pipeline_with_defaults\"></a>\n",
    "## Run Pipeline With Default Configuration\n",
    "\n",
    "Pipelines can be run by using the .call() method on the Pipeline class and passing in a data file.   Running a pipeline generally executes each successive step on the output from the previous step.  The end result is an output data model."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 39,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:44:12,779 - stpipe.Detector1Pipeline - INFO - Detector1Pipeline instance created.\n",
      "2019-09-09 18:44:12,782 - stpipe.Detector1Pipeline.group_scale - INFO - GroupScaleStep instance created.\n",
      "2019-09-09 18:44:12,784 - stpipe.Detector1Pipeline.dq_init - INFO - DQInitStep instance created.\n",
      "2019-09-09 18:44:12,785 - stpipe.Detector1Pipeline.saturation - INFO - SaturationStep instance created.\n",
      "2019-09-09 18:44:12,788 - stpipe.Detector1Pipeline.ipc - INFO - IPCStep instance created.\n",
      "2019-09-09 18:44:12,789 - stpipe.Detector1Pipeline.superbias - INFO - SuperBiasStep instance created.\n",
      "2019-09-09 18:44:12,792 - stpipe.Detector1Pipeline.refpix - INFO - RefPixStep instance created.\n",
      "2019-09-09 18:44:12,794 - stpipe.Detector1Pipeline.rscd - INFO - RSCD_Step instance created.\n",
      "2019-09-09 18:44:12,796 - stpipe.Detector1Pipeline.firstframe - INFO - FirstFrameStep instance created.\n",
      "2019-09-09 18:44:12,797 - stpipe.Detector1Pipeline.lastframe - INFO - LastFrameStep instance created.\n",
      "2019-09-09 18:44:12,800 - stpipe.Detector1Pipeline.linearity - INFO - LinearityStep instance created.\n",
      "2019-09-09 18:44:12,803 - stpipe.Detector1Pipeline.dark_current - INFO - DarkCurrentStep instance created.\n",
      "2019-09-09 18:44:12,805 - stpipe.Detector1Pipeline.persistence - INFO - PersistenceStep instance created.\n",
      "2019-09-09 18:44:12,808 - stpipe.Detector1Pipeline.jump - INFO - JumpStep instance created.\n",
      "2019-09-09 18:44:12,810 - stpipe.Detector1Pipeline.ramp_fit - INFO - RampFitStep instance created.\n",
      "2019-09-09 18:44:12,812 - stpipe.Detector1Pipeline.gain_scale - INFO - GainScaleStep instance created.\n",
      "2019-09-09 18:44:12,903 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline running with args ('test.fits',).\n",
      "2019-09-09 18:44:13,803 - stpipe.Detector1Pipeline - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/datamodels/util.py:165: NoTypeWarning: model_type not found. Opening test.fits as a RampModel\n",
      "  warnings.warn(errmsg, NoTypeWarning)\n",
      "\n",
      "2019-09-09 18:44:13,804 - stpipe.Detector1Pipeline - INFO - Prefetching reference files for dataset: 'test.fits' reftypes = ['dark', 'gain', 'ipc', 'linearity', 'mask', 'persat', 'readnoise', 'refpix', 'rscd', 'saturation', 'superbias', 'trapdensity', 'trappars']\n",
      "2019-09-09 18:44:14,951 - stpipe.Detector1Pipeline - INFO - Prefetch for DARK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits'.\n",
      "2019-09-09 18:44:14,952 - stpipe.Detector1Pipeline - INFO - Prefetch for GAIN reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits'.\n",
      "2019-09-09 18:44:14,953 - stpipe.Detector1Pipeline - INFO - Prefetch for IPC reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_ipc_0011.fits'.\n",
      "2019-09-09 18:44:14,954 - stpipe.Detector1Pipeline - INFO - Prefetch for LINEARITY reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits'.\n",
      "2019-09-09 18:44:14,955 - stpipe.Detector1Pipeline - INFO - Prefetch for MASK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits'.\n",
      "2019-09-09 18:44:14,956 - stpipe.Detector1Pipeline - INFO - Prefetch for PERSAT reference file is 'N/A'.\n",
      "2019-09-09 18:44:14,956 - stpipe.Detector1Pipeline - INFO - Prefetch for READNOISE reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits'.\n",
      "2019-09-09 18:44:14,957 - stpipe.Detector1Pipeline - INFO - Prefetch for REFPIX reference file is 'N/A'.\n",
      "2019-09-09 18:44:14,958 - stpipe.Detector1Pipeline - INFO - Prefetch for RSCD reference file is 'N/A'.\n",
      "2019-09-09 18:44:14,959 - stpipe.Detector1Pipeline - INFO - Prefetch for SATURATION reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits'.\n",
      "2019-09-09 18:44:14,960 - stpipe.Detector1Pipeline - INFO - Prefetch for SUPERBIAS reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits'.\n",
      "2019-09-09 18:44:14,961 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPDENSITY reference file is 'N/A'.\n",
      "2019-09-09 18:44:14,961 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPPARS reference file is 'N/A'.\n",
      "2019-09-09 18:44:14,963 - stpipe.Detector1Pipeline - INFO - Starting calwebb_detector1 ...\n",
      "2019-09-09 18:44:15,301 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:15,319 - stpipe.Detector1Pipeline.group_scale - INFO - NFRAMES and FRMDIVSR are equal; correction not needed\n",
      "2019-09-09 18:44:15,320 - stpipe.Detector1Pipeline.group_scale - INFO - Step will be skipped\n",
      "2019-09-09 18:44:15,322 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale done\n",
      "2019-09-09 18:44:15,387 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:15,407 - stpipe.Detector1Pipeline.dq_init - INFO - Using MASK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits\n",
      "2019-09-09 18:44:15,679 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init done\n",
      "2019-09-09 18:44:15,761 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:15,780 - stpipe.Detector1Pipeline.saturation - INFO - Using SATURATION reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits\n",
      "2019-09-09 18:44:16,058 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation done\n",
      "2019-09-09 18:44:16,124 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:16,149 - stpipe.Detector1Pipeline.ipc - INFO - Using IPC reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_ipc_0011.fits\n",
      "2019-09-09 18:44:16,650 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc done\n",
      "2019-09-09 18:44:16,714 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:16,734 - stpipe.Detector1Pipeline.superbias - INFO - Using SUPERBIAS reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits\n",
      "2019-09-09 18:44:16,945 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias done\n",
      "2019-09-09 18:44:17,013 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:17,027 - stpipe.Detector1Pipeline.refpix - INFO - use_side_ref_pixels = True\n",
      "2019-09-09 18:44:17,028 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_columns = True\n",
      "2019-09-09 18:44:17,029 - stpipe.Detector1Pipeline.refpix - INFO - side_smoothing_length = 11\n",
      "2019-09-09 18:44:17,031 - stpipe.Detector1Pipeline.refpix - INFO - side_gain = 1.000000\n",
      "2019-09-09 18:44:17,032 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_rows = True\n",
      "2019-09-09 18:44:18,338 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix done\n",
      "2019-09-09 18:44:18,407 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:18,426 - stpipe.Detector1Pipeline.linearity - INFO - Using Linearity reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits\n",
      "2019-09-09 18:44:18,577 - stpipe.Detector1Pipeline.linearity - WARNING - Keyword BAD_LIN_CORR does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:44:18,972 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity done\n",
      "2019-09-09 18:44:19,042 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:19,062 - stpipe.Detector1Pipeline.dark_current - INFO - Using DARK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits\n",
      "2019-09-09 18:44:19,393 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RTN does not correspond to an existing DQ mnemonic, so will be ignored\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:44:19,394 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RC_PIXEL does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:44:19,398 - stpipe.Detector1Pipeline.dark_current - INFO - Science data nints=1, ngroups=3, nframes=4, groupgap=0\n",
      "2019-09-09 18:44:19,399 - stpipe.Detector1Pipeline.dark_current - INFO - Dark data nints=1, ngroups=88, nframes=1, groupgap=0\n",
      "2019-09-09 18:44:21,182 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current done\n",
      "2019-09-09 18:44:21,244 - stpipe.Detector1Pipeline.jump - INFO - Step jump running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:21,261 - stpipe.Detector1Pipeline.jump - INFO - CR rejection threshold = 4 sigma\n",
      "2019-09-09 18:44:21,267 - stpipe.Detector1Pipeline.jump - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:44:21,325 - stpipe.Detector1Pipeline.jump - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:44:21,522 - stpipe.Detector1Pipeline.jump - INFO - Executing two-point difference method\n",
      "2019-09-09 18:44:21,681 - stpipe.Detector1Pipeline.jump - INFO -  working on integration 1\n",
      "2019-09-09 18:44:22,509 - stpipe.Detector1Pipeline.jump - INFO - From highest outlier Two point found 5938 pixels with at least one CR\n",
      "2019-09-09 18:44:22,646 - stpipe.Detector1Pipeline.jump - INFO - The execution time in seconds: 1.385386\n",
      "2019-09-09 18:44:22,654 - stpipe.Detector1Pipeline.jump - INFO - Step jump done\n",
      "2019-09-09 18:44:22,712 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:44:22,733 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:44:22,817 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:44:22,845 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using algorithm = ols\n",
      "2019-09-09 18:44:22,846 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using weighting = optimal\n",
      "2019-09-09 18:44:22,847 - stpipe.Detector1Pipeline.ramp_fit - INFO - Effective integration time per group: 42.94708\n",
      "2019-09-09 18:45:15,553 - stpipe.Detector1Pipeline.ramp_fit - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/ramp_fitting/ramp_fit.py:500: RuntimeWarning: invalid value encountered in multiply\n",
      "  var_p4[num_int,:,:,:] *= ( segs_4[num_int,:,:,:] > 0)\n",
      "\n",
      "2019-09-09 18:45:16,258 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of groups per integration: 3\n",
      "2019-09-09 18:45:16,259 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of integrations: 1\n",
      "2019-09-09 18:45:16,402 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit done\n",
      "2019-09-09 18:45:16,466 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale running with args (<ImageModel(2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:16,555 - stpipe.Detector1Pipeline.gain_scale - INFO - Rescaling by 1.0\n",
      "2019-09-09 18:45:16,567 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale done\n",
      "2019-09-09 18:45:16,568 - stpipe.Detector1Pipeline - INFO - ... ending calwebb_detector1\n",
      "2019-09-09 18:45:16,570 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline done\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<ImageModel(2048, 2048) from test.fits>"
      ]
     },
     "execution_count": 39,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "output_model = Detector1Pipeline.call(\"test.fits\")\n",
    "output_model"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id='pipeline_configs'></a>\n",
    "## About Configuration Files\n",
    "\n",
    "Configuration files are optional inputs for each step of the pipeline, as well as for the pipeline itself. These files list step-specific parameters, and can also be used to control which steps are run as part of the pipeline.\n",
    "\n",
    "You can get the full compliment of configuration files using the collect_pipeline_cfgs convenience function from the command line:\n",
    "\n",
    "    $ collect_pipeline_cfgs ./\n",
    "    \n",
    "The above can be executed in a notebook terminal window or by using the cell ! shell escape command.\n",
    "    \n",
    "This creates a copy of all configuration files, for all steps and all JWST Instruments. Note that default parameters in the config files are not necessarily optimized for any particular instrument.\n",
    "\n",
    "Each of these configuration files can be customized to control pipeline behavior. For example, the configuration file for the Level 1 detector pipeline is called calwebb_detector1.cfg and contains a list (not necessarily in order) of the steps run as part of the Level 1 pipeline.\n",
    "\n",
    "name = \"Detector1Pipeline\"\n",
    "class = \"jwst.pipeline.Detector1Pipeline\"\n",
    "save_calibrated_ramp = False\n",
    "\n",
    "    [steps]\n",
    "      [[group_scale]]\n",
    "        config_file = group_scale.cfg\n",
    "      [[dq_init]]\n",
    "        config_file = dq_init.cfg\n",
    "      [[saturation]]\n",
    "        config_file = saturation.cfg\n",
    "      [[ipc]]\n",
    "        skip = True\n",
    "      [[superbias]]\n",
    "        config_file = superbias.cfg\n",
    "      [[refpix]]\n",
    "        config_file = refpix.cfg\n",
    "      [[rscd]]\n",
    "        config_file = rscd.cfg\n",
    "      [[firstframe]]\n",
    "        config_file = firstframe.cfg\n",
    "      [[lastframe]]\n",
    "        config_file = lastframe.cfg\n",
    "      [[linearity]]\n",
    "        config_file = linearity.cfg\n",
    "      [[dark_current]]\n",
    "        config_file = dark_current.cfg\n",
    "      [[persistence]]\n",
    "        config_file = persistence.cfg\n",
    "      [[jump]]\n",
    "        config_file = jump.cfg\n",
    "      [[ramp_fit]]\n",
    "        config_file = ramp_fit.cfg\n",
    "      [[gain_scale]]\n",
    "        config_file = gain_scale.cfg\n",
    "\n",
    "In this example, the ipc step will be skipped (skip = True), and the output from the superbias step will be saved (save_results = True).\n",
    "\n",
    "Note that calwebb_detector1.cfg lists a configuration file for each pipeline step. You can customize a particular pipeline step by editing the parameters in its configuration file. \n",
    "\n",
    "For example, the persistence configuration file, shown below, contains details a flag *save_persistence* which when set to True will result in an additional output file with suffix _output_pers.\n",
    "\n",
    "name = \"persistence\"\n",
    "class = \"jwst.persistence.PersistenceStep\"\n",
    "\n",
    "input_trapsfilled = \"\"\n",
    "flag_pers_cutoff = 40.\n",
    "save_persistence = False\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Config File Setup\n",
    "\n",
    "To obtain all config files in a \"cfgs\" directory execute the following cell.   \n",
    "\n",
    "If you want to start over,  uncomment the \"rm -rf ./cfgs\" command.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 40,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "align_refs.cfg               guider_cds.cfg\r\n",
      "ami_analyze.cfg              hlsp.cfg\r\n",
      "ami_average.cfg              imprint.cfg\r\n",
      "ami_normalize.cfg            ipc.cfg\r\n",
      "assign_wcs.cfg               jump.cfg\r\n",
      "background.cfg               klip.cfg\r\n",
      "calwebb_ami3.cfg             lastframe.cfg\r\n",
      "calwebb_coron3.cfg           linear_pipeline.cfg\r\n",
      "calwebb_dark.cfg             linearity.cfg\r\n",
      "calwebb_detector1.cfg        master_background.cfg\r\n",
      "calwebb_detector1.cfg~       mrs_imatch.cfg\r\n",
      "calwebb_guider.cfg           outlier_detection.cfg\r\n",
      "calwebb_image2.cfg           outlier_detection_scaled.cfg\r\n",
      "calwebb_image3.cfg           outlier_detection_stack.cfg\r\n",
      "calwebb_nrslamp-spec2.cfg    outlier_detection_tso.cfg\r\n",
      "calwebb_spec2.cfg            pathloss.cfg\r\n",
      "calwebb_spec3.cfg            persistence.cfg\r\n",
      "calwebb_tso-image2.cfg       persistence.cfg~\r\n",
      "calwebb_tso-spec2.cfg        photom.cfg\r\n",
      "calwebb_tso1.cfg             ramp_fit.cfg\r\n",
      "calwebb_tso3.cfg             refpix.cfg\r\n",
      "calwebb_wfs-image2.cfg       resample.cfg\r\n",
      "calwebb_wfs-image3.cfg       resample_spec.cfg\r\n",
      "combine_1d.cfg               reset.cfg\r\n",
      "cube_build.cfg               rscd.cfg\r\n",
      "dark_current.cfg             saturation.cfg\r\n",
      "dq_init.cfg                  skymatch.cfg\r\n",
      "emission.cfg                 source_catalog.cfg\r\n",
      "engdblog.cfg                 srctype.cfg\r\n",
      "extract_1d.cfg               stack_refs.cfg\r\n",
      "extract_2d.cfg               straylight.cfg\r\n",
      "firstframe.cfg               subtract_images.cfg\r\n",
      "flat_field.cfg               superbias.cfg\r\n",
      "fringe.cfg                   tweakreg.cfg\r\n",
      "gain_scale.cfg               wfs_combine.cfg\r\n",
      "group_scale.cfg              white_light.cfg\r\n"
     ]
    }
   ],
   "source": [
    "# ! rm -rf ./cfgs\n",
    "! [ ! -d \"./cfgs\" ] && collect_pipeline_cfgs cfgs\n",
    "! ls ./cfgs"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id=\"pipeline_with_cfgs\"></a>\n",
    "## Run Pipeline with Configuration Files\n",
    "\n",
    "In your cfgs directory, edit the file calwebb_detector1.cfg and change this:\n",
    "\n",
    "      [[dark_current]]\n",
    "        config_file = dark_current.cfg\n",
    "        \n",
    "to this:\n",
    "\n",
    "      [[dark_current]]\n",
    "        config_file = dark_current.cfg\n",
    "        save_results = True\n",
    "\n",
    "This will command the DarkCurrentStep to save its output in a file which can be used to run the PersistenceStep in a standalone mode later.\n",
    "\n",
    "Then call the Detector1Pipeline with your modified config file."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 41,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:45:16,858 - stpipe.Detector1Pipeline - INFO - Detector1Pipeline instance created.\n",
      "2019-09-09 18:45:16,860 - stpipe.Detector1Pipeline.group_scale - INFO - GroupScaleStep instance created.\n",
      "2019-09-09 18:45:16,862 - stpipe.Detector1Pipeline.dq_init - INFO - DQInitStep instance created.\n",
      "2019-09-09 18:45:16,864 - stpipe.Detector1Pipeline.saturation - INFO - SaturationStep instance created.\n",
      "2019-09-09 18:45:16,866 - stpipe.Detector1Pipeline.ipc - INFO - IPCStep instance created.\n",
      "2019-09-09 18:45:16,868 - stpipe.Detector1Pipeline.superbias - INFO - SuperBiasStep instance created.\n",
      "2019-09-09 18:45:16,870 - stpipe.Detector1Pipeline.refpix - INFO - RefPixStep instance created.\n",
      "2019-09-09 18:45:16,872 - stpipe.Detector1Pipeline.rscd - INFO - RSCD_Step instance created.\n",
      "2019-09-09 18:45:16,874 - stpipe.Detector1Pipeline.firstframe - INFO - FirstFrameStep instance created.\n",
      "2019-09-09 18:45:16,877 - stpipe.Detector1Pipeline.lastframe - INFO - LastFrameStep instance created.\n",
      "2019-09-09 18:45:16,880 - stpipe.Detector1Pipeline.linearity - INFO - LinearityStep instance created.\n",
      "2019-09-09 18:45:16,882 - stpipe.Detector1Pipeline.dark_current - INFO - DarkCurrentStep instance created.\n",
      "2019-09-09 18:45:16,883 - stpipe.Detector1Pipeline.persistence - INFO - PersistenceStep instance created.\n",
      "2019-09-09 18:45:16,886 - stpipe.Detector1Pipeline.jump - INFO - JumpStep instance created.\n",
      "2019-09-09 18:45:16,888 - stpipe.Detector1Pipeline.ramp_fit - INFO - RampFitStep instance created.\n",
      "2019-09-09 18:45:16,890 - stpipe.Detector1Pipeline.gain_scale - INFO - GainScaleStep instance created.\n",
      "2019-09-09 18:45:16,971 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline running with args ('test.fits',).\n",
      "2019-09-09 18:45:17,269 - stpipe.Detector1Pipeline - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/datamodels/util.py:165: NoTypeWarning: model_type not found. Opening test.fits as a RampModel\n",
      "  warnings.warn(errmsg, NoTypeWarning)\n",
      "\n",
      "2019-09-09 18:45:17,271 - stpipe.Detector1Pipeline - INFO - Prefetching reference files for dataset: 'test.fits' reftypes = ['dark', 'gain', 'linearity', 'mask', 'persat', 'readnoise', 'refpix', 'rscd', 'saturation', 'superbias', 'trapdensity', 'trappars']\n",
      "2019-09-09 18:45:17,276 - stpipe.Detector1Pipeline - INFO - Prefetch for DARK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits'.\n",
      "2019-09-09 18:45:17,277 - stpipe.Detector1Pipeline - INFO - Prefetch for GAIN reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits'.\n",
      "2019-09-09 18:45:17,278 - stpipe.Detector1Pipeline - INFO - Prefetch for LINEARITY reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits'.\n",
      "2019-09-09 18:45:17,279 - stpipe.Detector1Pipeline - INFO - Prefetch for MASK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits'.\n",
      "2019-09-09 18:45:17,280 - stpipe.Detector1Pipeline - INFO - Prefetch for PERSAT reference file is 'N/A'.\n",
      "2019-09-09 18:45:17,281 - stpipe.Detector1Pipeline - INFO - Prefetch for READNOISE reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits'.\n",
      "2019-09-09 18:45:17,282 - stpipe.Detector1Pipeline - INFO - Prefetch for REFPIX reference file is 'N/A'.\n",
      "2019-09-09 18:45:17,283 - stpipe.Detector1Pipeline - INFO - Prefetch for RSCD reference file is 'N/A'.\n",
      "2019-09-09 18:45:17,284 - stpipe.Detector1Pipeline - INFO - Prefetch for SATURATION reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits'.\n",
      "2019-09-09 18:45:17,285 - stpipe.Detector1Pipeline - INFO - Prefetch for SUPERBIAS reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits'.\n",
      "2019-09-09 18:45:17,286 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPDENSITY reference file is 'N/A'.\n",
      "2019-09-09 18:45:17,286 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPPARS reference file is 'N/A'.\n",
      "2019-09-09 18:45:17,287 - stpipe.Detector1Pipeline - INFO - Starting calwebb_detector1 ...\n",
      "2019-09-09 18:45:17,622 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:17,637 - stpipe.Detector1Pipeline.group_scale - INFO - NFRAMES and FRMDIVSR are equal; correction not needed\n",
      "2019-09-09 18:45:17,638 - stpipe.Detector1Pipeline.group_scale - INFO - Step will be skipped\n",
      "2019-09-09 18:45:17,640 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale done\n",
      "2019-09-09 18:45:17,701 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:17,722 - stpipe.Detector1Pipeline.dq_init - INFO - Using MASK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits\n",
      "2019-09-09 18:45:17,980 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init done\n",
      "2019-09-09 18:45:18,046 - stpipe.Detector1Pipeline - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/stpipe/step.py:341: ResourceWarning: unclosed file <_io.FileIO name='test.fits' mode='rb' closefd=True>\n",
      "  gc.collect()\n",
      "\n",
      "2019-09-09 18:45:18,070 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:18,087 - stpipe.Detector1Pipeline.saturation - INFO - Using SATURATION reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits\n",
      "2019-09-09 18:45:18,363 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation done\n",
      "2019-09-09 18:45:18,426 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:18,427 - stpipe.Detector1Pipeline.ipc - INFO - Step skipped.\n",
      "2019-09-09 18:45:18,430 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc done\n",
      "2019-09-09 18:45:18,494 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:18,514 - stpipe.Detector1Pipeline.superbias - INFO - Using SUPERBIAS reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits\n",
      "2019-09-09 18:45:18,716 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias done\n",
      "2019-09-09 18:45:18,783 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:18,798 - stpipe.Detector1Pipeline.refpix - INFO - use_side_ref_pixels = True\n",
      "2019-09-09 18:45:18,798 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_columns = True\n",
      "2019-09-09 18:45:18,799 - stpipe.Detector1Pipeline.refpix - INFO - side_smoothing_length = 11\n",
      "2019-09-09 18:45:18,800 - stpipe.Detector1Pipeline.refpix - INFO - side_gain = 1.000000\n",
      "2019-09-09 18:45:18,801 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_rows = True\n",
      "2019-09-09 18:45:20,182 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix done\n",
      "2019-09-09 18:45:20,244 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:20,280 - stpipe.Detector1Pipeline.linearity - INFO - Using Linearity reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits\n",
      "2019-09-09 18:45:20,406 - stpipe.Detector1Pipeline.linearity - WARNING - Keyword BAD_LIN_CORR does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:45:20,794 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity done\n",
      "2019-09-09 18:45:20,863 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:45:20,885 - stpipe.Detector1Pipeline.dark_current - INFO - Using DARK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits\n",
      "2019-09-09 18:45:21,190 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RTN does not correspond to an existing DQ mnemonic, so will be ignored\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:45:21,191 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RC_PIXEL does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:45:21,198 - stpipe.Detector1Pipeline.dark_current - INFO - Science data nints=1, ngroups=3, nframes=4, groupgap=0\n",
      "2019-09-09 18:45:21,198 - stpipe.Detector1Pipeline.dark_current - INFO - Dark data nints=1, ngroups=88, nframes=1, groupgap=0\n",
      "2019-09-09 18:45:23,593 - stpipe.Detector1Pipeline.dark_current - INFO - Saved model in test_dark_current.fits\n",
      "2019-09-09 18:45:23,594 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current done\n",
      "2019-09-09 18:45:23,655 - stpipe.Detector1Pipeline.jump - INFO - Step jump running with args (<RampModel(1, 3, 2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:45:23,689 - stpipe.Detector1Pipeline.jump - INFO - CR rejection threshold = 4 sigma\n",
      "2019-09-09 18:45:23,694 - stpipe.Detector1Pipeline.jump - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:45:23,724 - stpipe.Detector1Pipeline.jump - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:45:23,925 - stpipe.Detector1Pipeline.jump - INFO - Executing two-point difference method\n",
      "2019-09-09 18:45:24,092 - stpipe.Detector1Pipeline.jump - INFO -  working on integration 1\n",
      "2019-09-09 18:45:24,946 - stpipe.Detector1Pipeline.jump - INFO - From highest outlier Two point found 5120 pixels with at least one CR\n",
      "2019-09-09 18:45:25,080 - stpipe.Detector1Pipeline.jump - INFO - The execution time in seconds: 1.391092\n",
      "2019-09-09 18:45:25,089 - stpipe.Detector1Pipeline.jump - INFO - Step jump done\n",
      "2019-09-09 18:45:25,149 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit running with args (<RampModel(1, 3, 2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:45:25,179 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:45:25,246 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:45:25,264 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using algorithm = ols\n",
      "2019-09-09 18:45:25,264 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using weighting = optimal\n",
      "2019-09-09 18:45:25,265 - stpipe.Detector1Pipeline.ramp_fit - INFO - Effective integration time per group: 42.94708\n",
      "2019-09-09 18:46:08,754 - stpipe.Detector1Pipeline.ramp_fit - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/ramp_fitting/ramp_fit.py:500: RuntimeWarning: invalid value encountered in multiply\n",
      "  var_p4[num_int,:,:,:] *= ( segs_4[num_int,:,:,:] > 0)\n",
      "\n",
      "2019-09-09 18:46:09,425 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of groups per integration: 3\n",
      "2019-09-09 18:46:09,427 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of integrations: 1\n",
      "2019-09-09 18:46:09,529 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit done\n",
      "2019-09-09 18:46:09,583 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale running with args (<ImageModel(2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:46:09,646 - stpipe.Detector1Pipeline.gain_scale - INFO - Rescaling by 1.0\n",
      "2019-09-09 18:46:09,656 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale done\n",
      "2019-09-09 18:46:09,657 - stpipe.Detector1Pipeline - INFO - ... ending calwebb_detector1\n",
      "2019-09-09 18:46:09,659 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline done\n"
     ]
    }
   ],
   "source": [
    "output_model = Detector1Pipeline.call(\"test.fits\", config_file='cfgs/calwebb_detector1.cfg')\n",
    "output_model"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 43,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "test_dark_current.fits\r\n"
     ]
    }
   ],
   "source": [
    "!ls test_dark_current.fits   # you should see a file since you set save_results=True"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 53,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 864x864 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_image(output_model.data)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 54,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 864x864 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_image(\"test_dark_current.fits\", index=(0,0))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id=\"pipeline_no_configs\"></a>\n",
    "## Run Pipeline with Parameters Set Programmatically\n",
    "\n",
    "You can also run the pipeline without relying on configuration files by setting parameters programmatically, and relying on the defaults in the pipeline."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 56,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:55:38,357 - stpipe.Detector1Pipeline - INFO - Detector1Pipeline instance created.\n",
      "2019-09-09 18:55:38,359 - stpipe.Detector1Pipeline.group_scale - INFO - GroupScaleStep instance created.\n",
      "2019-09-09 18:55:38,361 - stpipe.Detector1Pipeline.dq_init - INFO - DQInitStep instance created.\n",
      "2019-09-09 18:55:38,363 - stpipe.Detector1Pipeline.saturation - INFO - SaturationStep instance created.\n",
      "2019-09-09 18:55:38,365 - stpipe.Detector1Pipeline.ipc - INFO - IPCStep instance created.\n",
      "2019-09-09 18:55:38,367 - stpipe.Detector1Pipeline.superbias - INFO - SuperBiasStep instance created.\n",
      "2019-09-09 18:55:38,369 - stpipe.Detector1Pipeline.refpix - INFO - RefPixStep instance created.\n",
      "2019-09-09 18:55:38,371 - stpipe.Detector1Pipeline.rscd - INFO - RSCD_Step instance created.\n",
      "2019-09-09 18:55:38,373 - stpipe.Detector1Pipeline.firstframe - INFO - FirstFrameStep instance created.\n",
      "2019-09-09 18:55:38,375 - stpipe.Detector1Pipeline.lastframe - INFO - LastFrameStep instance created.\n",
      "2019-09-09 18:55:38,377 - stpipe.Detector1Pipeline.linearity - INFO - LinearityStep instance created.\n",
      "2019-09-09 18:55:38,379 - stpipe.Detector1Pipeline.dark_current - INFO - DarkCurrentStep instance created.\n",
      "2019-09-09 18:55:38,381 - stpipe.Detector1Pipeline.persistence - INFO - PersistenceStep instance created.\n",
      "2019-09-09 18:55:38,383 - stpipe.Detector1Pipeline.jump - INFO - JumpStep instance created.\n",
      "2019-09-09 18:55:38,385 - stpipe.Detector1Pipeline.ramp_fit - INFO - RampFitStep instance created.\n",
      "2019-09-09 18:55:38,387 - stpipe.Detector1Pipeline.gain_scale - INFO - GainScaleStep instance created.\n",
      "2019-09-09 18:55:38,529 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline running with args ('test.fits',).\n",
      "2019-09-09 18:55:38,822 - stpipe.Detector1Pipeline - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/datamodels/util.py:165: NoTypeWarning: model_type not found. Opening test.fits as a RampModel\n",
      "  warnings.warn(errmsg, NoTypeWarning)\n",
      "\n",
      "2019-09-09 18:55:38,823 - stpipe.Detector1Pipeline - INFO - Prefetching reference files for dataset: 'test.fits' reftypes = ['dark', 'gain', 'ipc', 'linearity', 'mask', 'persat', 'readnoise', 'refpix', 'rscd', 'saturation', 'superbias', 'trapdensity', 'trappars']\n",
      "2019-09-09 18:55:38,829 - stpipe.Detector1Pipeline - INFO - Prefetch for DARK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits'.\n",
      "2019-09-09 18:55:38,829 - stpipe.Detector1Pipeline - INFO - Prefetch for GAIN reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits'.\n",
      "2019-09-09 18:55:38,830 - stpipe.Detector1Pipeline - INFO - Prefetch for IPC reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_ipc_0011.fits'.\n",
      "2019-09-09 18:55:38,831 - stpipe.Detector1Pipeline - INFO - Prefetch for LINEARITY reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits'.\n",
      "2019-09-09 18:55:38,832 - stpipe.Detector1Pipeline - INFO - Prefetch for MASK reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits'.\n",
      "2019-09-09 18:55:38,832 - stpipe.Detector1Pipeline - INFO - Prefetch for PERSAT reference file is 'N/A'.\n",
      "2019-09-09 18:55:38,833 - stpipe.Detector1Pipeline - INFO - Prefetch for READNOISE reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits'.\n",
      "2019-09-09 18:55:38,834 - stpipe.Detector1Pipeline - INFO - Prefetch for REFPIX reference file is 'N/A'.\n",
      "2019-09-09 18:55:38,835 - stpipe.Detector1Pipeline - INFO - Prefetch for RSCD reference file is 'N/A'.\n",
      "2019-09-09 18:55:38,835 - stpipe.Detector1Pipeline - INFO - Prefetch for SATURATION reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits'.\n",
      "2019-09-09 18:55:38,836 - stpipe.Detector1Pipeline - INFO - Prefetch for SUPERBIAS reference file is '/Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits'.\n",
      "2019-09-09 18:55:38,837 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPDENSITY reference file is 'N/A'.\n",
      "2019-09-09 18:55:38,838 - stpipe.Detector1Pipeline - INFO - Prefetch for TRAPPARS reference file is 'N/A'.\n",
      "2019-09-09 18:55:38,839 - stpipe.Detector1Pipeline - INFO - Starting calwebb_detector1 ...\n",
      "2019-09-09 18:55:39,186 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:39,203 - stpipe.Detector1Pipeline.group_scale - INFO - NFRAMES and FRMDIVSR are equal; correction not needed\n",
      "2019-09-09 18:55:39,204 - stpipe.Detector1Pipeline.group_scale - INFO - Step will be skipped\n",
      "2019-09-09 18:55:39,206 - stpipe.Detector1Pipeline.group_scale - INFO - Step group_scale done\n",
      "2019-09-09 18:55:39,275 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:39,298 - stpipe.Detector1Pipeline.dq_init - INFO - Using MASK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_mask_0024.fits\n",
      "2019-09-09 18:55:39,566 - stpipe.Detector1Pipeline.dq_init - INFO - Step dq_init done\n",
      "2019-09-09 18:55:39,630 - stpipe.Detector1Pipeline - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/stpipe/step.py:341: ResourceWarning: unclosed file <_io.FileIO name='test.fits' mode='rb' closefd=True>\n",
      "  gc.collect()\n",
      "\n",
      "2019-09-09 18:55:39,657 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:39,674 - stpipe.Detector1Pipeline.saturation - INFO - Using SATURATION reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_saturation_0027.fits\n",
      "2019-09-09 18:55:39,941 - stpipe.Detector1Pipeline.saturation - INFO - Step saturation done\n",
      "2019-09-09 18:55:40,009 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:40,035 - stpipe.Detector1Pipeline.ipc - INFO - Using IPC reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_ipc_0011.fits\n",
      "2019-09-09 18:55:40,523 - stpipe.Detector1Pipeline.ipc - INFO - Step ipc done\n",
      "2019-09-09 18:55:40,588 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:40,608 - stpipe.Detector1Pipeline.superbias - INFO - Using SUPERBIAS reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_superbias_0189.fits\n",
      "2019-09-09 18:55:40,822 - stpipe.Detector1Pipeline.superbias - INFO - Step superbias done\n",
      "2019-09-09 18:55:40,888 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:40,904 - stpipe.Detector1Pipeline.refpix - INFO - use_side_ref_pixels = True\n",
      "2019-09-09 18:55:40,905 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_columns = True\n",
      "2019-09-09 18:55:40,906 - stpipe.Detector1Pipeline.refpix - INFO - side_smoothing_length = 11\n",
      "2019-09-09 18:55:40,907 - stpipe.Detector1Pipeline.refpix - INFO - side_gain = 1.000000\n",
      "2019-09-09 18:55:40,908 - stpipe.Detector1Pipeline.refpix - INFO - odd_even_rows = True\n",
      "2019-09-09 18:55:42,268 - stpipe.Detector1Pipeline.refpix - INFO - Step refpix done\n",
      "2019-09-09 18:55:42,340 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:42,358 - stpipe.Detector1Pipeline.linearity - INFO - Using Linearity reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_linearity_0022.fits\n",
      "2019-09-09 18:55:42,499 - stpipe.Detector1Pipeline.linearity - WARNING - Keyword BAD_LIN_CORR does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:55:42,876 - stpipe.Detector1Pipeline.linearity - INFO - Step linearity done\n",
      "2019-09-09 18:55:42,943 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current running with args (<RampModel(1, 3, 2048, 2048) from test.fits>,).\n",
      "2019-09-09 18:55:42,963 - stpipe.Detector1Pipeline.dark_current - INFO - Using DARK reference file /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_dark_0151.fits\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:55:43,276 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RTN does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:55:43,277 - stpipe.Detector1Pipeline.dark_current - WARNING - Keyword RC_PIXEL does not correspond to an existing DQ mnemonic, so will be ignored\n",
      "2019-09-09 18:55:43,284 - stpipe.Detector1Pipeline.dark_current - INFO - Science data nints=1, ngroups=3, nframes=4, groupgap=0\n",
      "2019-09-09 18:55:43,285 - stpipe.Detector1Pipeline.dark_current - INFO - Dark data nints=1, ngroups=88, nframes=1, groupgap=0\n",
      "2019-09-09 18:55:45,593 - stpipe.Detector1Pipeline.dark_current - INFO - Saved model in test_dark_current.fits\n",
      "2019-09-09 18:55:45,594 - stpipe.Detector1Pipeline.dark_current - INFO - Step dark_current done\n",
      "2019-09-09 18:55:45,656 - stpipe.Detector1Pipeline.jump - INFO - Step jump running with args (<RampModel(1, 3, 2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:55:45,688 - stpipe.Detector1Pipeline.jump - INFO - CR rejection threshold = 4 sigma\n",
      "2019-09-09 18:55:45,692 - stpipe.Detector1Pipeline.jump - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:55:45,727 - stpipe.Detector1Pipeline.jump - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:55:45,919 - stpipe.Detector1Pipeline.jump - INFO - Executing two-point difference method\n",
      "2019-09-09 18:55:46,081 - stpipe.Detector1Pipeline.jump - INFO -  working on integration 1\n",
      "2019-09-09 18:55:46,903 - stpipe.Detector1Pipeline.jump - INFO - From highest outlier Two point found 5938 pixels with at least one CR\n",
      "2019-09-09 18:55:47,041 - stpipe.Detector1Pipeline.jump - INFO - The execution time in seconds: 1.352160\n",
      "2019-09-09 18:55:47,048 - stpipe.Detector1Pipeline.jump - INFO - Step jump done\n",
      "2019-09-09 18:55:47,114 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit running with args (<RampModel(1, 3, 2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:55:47,143 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using READNOISE reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_readnoise_0019.fits\n",
      "2019-09-09 18:55:47,219 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using GAIN reference file: /Users/jmiller/crds_cache_ops/references/jwst/nirspec/jwst_nirspec_gain_0023.fits\n",
      "2019-09-09 18:55:47,247 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using algorithm = ols\n",
      "2019-09-09 18:55:47,248 - stpipe.Detector1Pipeline.ramp_fit - INFO - Using weighting = optimal\n",
      "2019-09-09 18:55:47,249 - stpipe.Detector1Pipeline.ramp_fit - INFO - Effective integration time per group: 42.94708\n",
      "2019-09-09 18:56:38,408 - stpipe.Detector1Pipeline.ramp_fit - WARNING - /Users/jmiller/miniconda3/envs/crds-env/lib/python3.7/site-packages/jwst/ramp_fitting/ramp_fit.py:500: RuntimeWarning: invalid value encountered in multiply\n",
      "  var_p4[num_int,:,:,:] *= ( segs_4[num_int,:,:,:] > 0)\n",
      "\n",
      "2019-09-09 18:56:39,087 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of groups per integration: 3\n",
      "2019-09-09 18:56:39,088 - stpipe.Detector1Pipeline.ramp_fit - INFO - Number of integrations: 1\n",
      "2019-09-09 18:56:39,208 - stpipe.Detector1Pipeline.ramp_fit - INFO - Step ramp_fit done\n",
      "2019-09-09 18:56:39,276 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale running with args (<ImageModel(2048, 2048) from test_dark_current.fits>,).\n",
      "2019-09-09 18:56:39,363 - stpipe.Detector1Pipeline.gain_scale - INFO - Rescaling by 1.0\n",
      "2019-09-09 18:56:39,374 - stpipe.Detector1Pipeline.gain_scale - INFO - Step gain_scale done\n",
      "2019-09-09 18:56:39,375 - stpipe.Detector1Pipeline - INFO - ... ending calwebb_detector1\n",
      "2019-09-09 18:56:39,377 - stpipe.Detector1Pipeline - INFO - Step Detector1Pipeline done\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<ImageModel(2048, 2048) from test_dark_current.fits>"
      ]
     },
     "execution_count": 56,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "! rm -f test_dark_current.fits   # remove the dark current output,  we'll recreate it another way...\n",
    "\n",
    "det1p = Detector1Pipeline()\n",
    "det1p.dark_current.save_results = True\n",
    "det1p.run(\"test.fits\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Setting save_results=True results in the output of the test_dark_current.fits file which can be passed into the PersistenceStep when run in isolation, next."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 57,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "test_dark_current.fits\r\n"
     ]
    }
   ],
   "source": [
    "! ls test_dark_current.fits   # make sure save_results worked."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "***\n",
    "<a id=\"steps_with_config_files\"></a>\n",
    "## Run Individual Steps with Configuration Files\n",
    "\n",
    "You can also change parameter values in the .cfg files for individual Steps.\n",
    "\n",
    "Edit the cfgs/persistence.cfg file and change:\n",
    "\n",
    "```\n",
    "   input_trapsfilled = \"\"\n",
    "   flag_pers_cutoff = 40.\n",
    "   save_persistence = False\n",
    "```\n",
    "\n",
    "to:\n",
    "\n",
    "```\n",
    "   input_trapsfilled = \"\"\n",
    "   flag_pers_cutoff = 40.\n",
    "   save_persistence = True\n",
    "```\n",
    "\n",
    "This will cause PersistenceStep to output a third output file with suffix “_output_pers”. \n",
    "\n",
    "See https://jwst-pipeline.readthedocs.io/en/latest/jwst/persistence/description.html\n",
    "for more information on the persistence step."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 58,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "2019-09-09 18:57:11,738 - stpipe.persistence - INFO - PersistenceStep instance created.\n",
      "2019-09-09 18:57:11,831 - stpipe.persistence - INFO - Step persistence running with args ('test_dark_current.fits',).\n",
      "2019-09-09 18:57:12,445 - stpipe.persistence - WARNING - Missing reference file types:  PERSAT TRAPDENSITY TRAPPARS\n",
      "2019-09-09 18:57:12,449 - stpipe.persistence - INFO - Step persistence done\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "<RampModel(1, 3, 2048, 2048) from test_dark_current.fits>"
      ]
     },
     "execution_count": 58,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "os.environ[\"\"]\n",
    "\n",
    "from jwst.persistence.persistence_step import PersistenceStep\n",
    "\n",
    "output_model = PersistenceStep.call(\"test_dark_current.fits\", config_file=\"cfgs/persistence.cfg\")\n",
    "\n",
    "output_model"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 62,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "iVBORw0KGgoAAAANSUhEUgAAAp8AAAKrCAYAAAC3A+azAAAABHNCSVQICAgIfAhkiAAAAAlwSFlzAAALEgAACxIB0t1+/AAAADh0RVh0U29mdHdhcmUAbWF0cGxvdGxpYiB2ZXJzaW9uMy4xLjEsIGh0dHA6Ly9tYXRwbG90bGliLm9yZy8QZhcZAAAgAElEQVR4nO3df7Bmd10n+PdnQggjiISJE2MSDLrB3WA5AbIxo6MVhxFCyiI44zjJbkFgWFvHsCu1bk2Bs2UsLKpcFa1l0VBRUiFbGED5ldmKYmTcQas2SBIzIT9kaBA23RWShTgEByeS7s/+cU/LQ6e7b3ffPud7+z6vV9WpPs/3nOc833Puc+/z6c/n+Z5vdXcAAGAJf2d0BwAAWB+CTwAAFiP4BABgMYJPAAAWI/gEAGAxgk8AABYj+AQA2Iaq6oaqeqSq7j2KfX+gqu6qqieq6kcP2vacqvqDqnqgqu6vqvPm6vPREHwCAGxPNya57Cj3/X+TvDrJbx9i201Jfrm7/5skFyd55ER07ngJPgEAtqHu/miSR1fbquo7qur3q+rOqvrjqvqvp30/2933JNl/0P4XJHlKd9827fdX3f2VhU7hkASfAAAnj+uT/I/d/aIk/0uS39hk/+cl+U9V9f6q+rOq+uWqOmX2Xh7BU0a+OAAAR6eqnpHke5P8TlUdaD5tk6c9Jcn3J3lBNkrz78lGef4d8/Ryc4JPAICTw99J8p+6+8JjeM6eJHd392eSpKo+mOSSDAw+ld0BAE4C3f1Ykr+oqn+eJLXhH2zytI8neVZVffP0+B8nuX/Gbm6qunvk6wMAcAhVdXOSS5OckeThJNcm+XdJrktyVpJTk7y7u99UVf9tkg8kOT3Jf0ny+e5+/nScH0ryliSV5M4ku7r7b5Y9m68RfAIAsBhldwAAFiP4BABgMUa7AwAM8NIffHp/8dF9w17/znse/3B3H+0MSieM4BMAYIAvProvf/rh5wx7/VPO+tQZI15X2R0AgMXIfAIADNBJ9n/9VOxrQeYTAIDFyHwCAAzR2dcynwAAMBvBJwAAi1F2BwAYYGPA0fpNcy7zCQDAYgSfAAAsRtkdAGAQ9/kEAIAZyXwCAAzQ6exrA44AAGA2gk8AABaj7A4AMIj7fAIAwIxkPgEABugk+2Q+AQBgPoJPAAAWo+wOADCIAUcAADAjmU8AgAE6McMRAADMSfAJAMBilN0BAAbZP7oDA8h8AgCwGJlPAIABOm2GIwAAmJPgEwCAxSi7AwCM0Mm+9au6y3wCALAcwScAAItRdgcAGKDjPp8AADArmU8AgCEq+1KjO7E4mU8AABYj+AQAYDHK7gAAA3SS/e7zCQAA85H5BAAYxIAjAACYkeATAIDFKLsDAAzQUXYHAIAkSVWdW1V/VFX3V9V9VfXTh9jn0qr6UlXdPS0/t9lxZT4BAAbZ39s68/lEkp/p7ruq6huT3FlVt3X3/Qft98fd/cNHe1CZTwAAnqS7H+ruu6b1Lyd5IMnZWz2u4BMAYD2dUVV3rCy7DrdjVZ2X5AVJPnaIzf+wqv5DVf1eVT1/sxdVdgcAGGAbDDj6QndftNlOVfWMJO9L8vrufuygzXcl+bbu/ququjzJB5Ocf6TjyXwCAHBIVXVqNgLPd3X3+w/e3t2PdfdfTeu3Jjm1qs440jFlPgEABuhU9m3jPGBVVZJ3JHmgu3/1MPt8S5KHu7ur6uJsJDa/eKTjCj4BADiU70vyyiSfqKq7p7afTfKcJOnutyf50ST/qqqeSPLXSa7s7j7SQQWfAAA8SXf/SXLkL6V299uSvO1Yjiv4BAAYZJvf53MW2/eLBgAA7DgynwAAA2yDWy0NIfMJAMBitn3m86l1Wj8tTx/dDQBgh/gv+c/5m358/VKO28S2Dz6flqfne+rFo7sBAOwQH+uPjO7CpLKv168IvX5nDADAMIJPAAAWs+3L7gAAO1En2b+GecD1O2MAAIaR+QQAGMR9PgEAYEaCTwAAFrNp8FlV51bVH1XV/VV1X1X99NT+7Kq6rao+Nf17+tReVfXWqtpdVfdU1QtXjnX1tP+nqurq+U4LAGB76964z+eoZZSjeeUnkvxMd1+Q5JIk11TVBUnekOQj3X1+ko9Mj5PkZUnOn5ZdSa5LNoLVJNcm+Z4kFye59kDACgDAetg0+Ozuh7r7rmn9y0keSHJ2kiuSvHPa7Z1JXjGtX5Hkpt5we5JnVdVZSV6a5LbufrS7/zLJbUkuO6FnAwBwEtmfGraMckw516o6L8kLknwsyZnd/dC06fNJzpzWz07y4MrT9kxth2sHAGBNHHXwWVXPSPK+JK/v7sdWt3V3Z+NeqSdEVe2qqjuq6o6v5vETdVgAAAY7qvt8VtWp2Qg839Xd75+aH66qs7r7oams/sjUvjfJuStPP2dq25vk0oPa/+9DvV53X5/k+iR5Zj37hAW1AADbRSfZt4Y3Hjqa0e6V5B1JHujuX13ZdEuSAyPWr07yoZX2V02j3i9J8qWpPP/hJC+pqtOngUYvmdoAAFgTR5P5/L4kr0zyiaq6e2r72SS/mOS9VfXaJJ9L8mPTtluTXJ5kd5KvJHlNknT3o1X1C0k+Pu33pu5+9IScBQDASaeG3vJolE2Dz+7+k+SwQ6JefIj9O8k1hznWDUluOJYOAgCwc6xfuA0AwDBHNeAIAIATq5PsX8M84PqdMQAAw8h8AgAMsq/HzTQ0iswnAACLEXwCALAYZXcAgAE6ZYYjAACYk8wnAMAg+9dwhqP1O2MAAIYRfAIAsBhldwCAATox4AgAAOYk+AQAYDHK7gAAA3TK9JoAADAnmU8AgEH2r2EecP3OGACAYQSfAAAsRtkdAGCA7mSf6TUBAGA+Mp8AAENU9setlgAAYDaCTwAAFqPsDgAwQMeAIwAAmJXMJwDAIPvWMA+4fmcMAMAwgk8AABaj7A4AMECnsr/d5xMAAGYj8wkAMIgBRwAAMCPBJwAAi1F2BwAYoJPsN8MRAADMR/AJAMBilN0BAIao7Iv7fAIAwGxkPgEABjDgCAAAZib4BABgMcruAACDGHAEAAAzkvkEABiguww4AgCAOQk+AQBYjLI7AMAg+5TdAQBgPjKfAAADdJL9brUEAADzEXwCALAYZXcAgCHKgCMAAJiTzCcAwACdZH8bcAQAALPZNPisqhuq6pGqunel7T1Vdfe0fLaq7p7az6uqv17Z9vaV57yoqj5RVbur6q1VtX6hPgDAmjuasvuNSd6W5KYDDd39Lw6sV9VbknxpZf9Pd/eFhzjOdUl+PMnHktya5LIkv3fsXQYA2Bn2rWERetMz7u6PJnn0UNum7OWPJbn5SMeoqrOSPLO7b+/uzkYg+4pj7y4AACezrQ44+v4kD3f3p1banltVf5bksST/a3f/cZKzk+xZ2WfP1AYAsJY6tZYDjrYafF6Vr896PpTkOd39xap6UZIPVtXzj/WgVbUrya4keVq+YYtdBABguzju4LOqnpLknyZ50YG27n48yePT+p1V9ekkz0uyN8k5K08/Z2o7pO6+Psn1SfLMenYfbx8BANhetpL5/CdJ/ry7/7acXlXfnOTR7t5XVd+e5Pwkn+nuR6vqsaq6JBsDjl6V5P/YSscBAE52+w04erKqujnJ/5PkO6tqT1W9dtp0ZZ480OgHktwz3Xrpd5P8ZHcfGKz0U0l+K8nuJJ+Oke4AAGtn08xnd191mPZXH6LtfUned5j970jyXcfYPwAAdhDTawIADNCd7FvD0e7r90UDAACGkfkEABhkHe/zKfMJAMBiBJ8AACxG2R0AYICN6TXXLw+4fmcMAMAwMp8AAIPsiwFHAAAwG8EnAACLUXYHABig4z6fAAAwK5lPAIAh3GoJAABmJfgEAGAxyu4AAIPsd59PAACYj8wnAMAA3ck+t1oCAID5CD4BAFiMsjsAwCDu8wkAADOS+QQAGKBT5nYHAIA5CT4BAHiSqjq3qv6oqu6vqvuq6qcPsU9V1VurandV3VNVL9zsuMruAACDbPMZjp5I8jPdfVdVfWOSO6vqtu6+f2WflyU5f1q+J8l107+HJfMJAMCTdPdD3X3XtP7lJA8kOfug3a5IclNvuD3Js6rqrCMdV/AJAMARVdV5SV6Q5GMHbTo7yYMrj/fkyQHq11F2BwAYoJPRo93PqKo7Vh5f393XH7xTVT0jyfuSvL67H9vqiwo+AQDW0xe6+6Ij7VBVp2Yj8HxXd7//ELvsTXLuyuNzprbDEnwCAAyynWc4qqpK8o4kD3T3rx5mt1uSvK6q3p2NgUZf6u6HjnRcwScAAIfyfUlemeQTVXX31PazSZ6TJN399iS3Jrk8ye4kX0nyms0OKvgEAOBJuvtPkiPfC6q7O8k1x3JcwScAwAhtek0AAJiVzCcAwACdbT/D0SxkPgEAWIzgEwCAxSi7AwAMYsARAADMSOYTAGCAbTC3+xAynwAALEbwCQDAYpTdAQAGUXYHAIAZyXwCAAzQMbc7AADMSvAJAMBilN0BAAbZH2V3AACYjeATAIDFKLsDAIzQ7vMJAACzkvkEABigI/MJAACzEnwCALAYZXcAgEGU3Q+hqm6oqkeq6t6Vtp+vqr1Vdfe0XL6y7Y1VtbuqPllVL11pv2xq211VbzjxpwIAwHZ3NJnPG5O8LclNB7X/Wnf/ympDVV2Q5Mokz0/yrUn+sKqeN23+9SQ/lGRPko9X1S3dff8W+g4AcNLq1FpmPjcNPrv7o1V13lEe74ok7+7ux5P8RVXtTnLxtG13d38mSarq3dO+gk8AgDWylQFHr6uqe6ay/OlT29lJHlzZZ8/Udrj2Q6qqXVV1R1Xd8dU8voUuAgCwnRxv8Hldku9IcmGSh5K85YT1KEl3X9/dF3X3RafmtBN5aACAbaO7hi2jHNdo9+5++MB6Vf1mkv9rerg3ybkru54zteUI7QAArInjynxW1VkrD38kyYGR8LckubKqTquq5yY5P8mfJvl4kvOr6rlV9dRsDEq65fi7DQBw8tufGraMsmnms6puTnJpkjOqak+Sa5NcWlUXZmNmqM8m+Ykk6e77quq92RhI9ESSa7p733Sc1yX5cJJTktzQ3fed8LMBAGBbO5rR7lcdovkdR9j/zUnefIj2W5Pceky9AwBgRzHDEQDAAN1mOAIAgFnJfAIADDLylkejyHwCALAYwScAAItRdgcAGKIMOAIAgDnJfAIADGLAEQAAzEjwCQDAYpTdAQAG6JjhCAAAZiX4BABgMcruAAAjdNI9uhPLk/kEAGAxMp8AAIPsjwFHAAAwG8EnAACLUXYHABigY3pNAACYlcwnAMAQZYYjAACYk+ATAIDFKLsDAAxihiMAAJiRzCcAwCButQQAADMSfAIAsBhldwCAAbqV3QEAYFYynwAAg5jhCAAAZiT4BABgMcruAACDmOEIAABmJPMJADCIWy0BAMCMBJ8AACxG2R0AYIBOKbsDAMCcBJ8AACxG2R0AYJA1vM2nzCcAAMuR+QQAGKHd5xMAAGYl+AQAYDHK7gAAo6zhiCOZTwAAFiPzCQAwiAFHAAAwI8EnAACLUXYHABikDTgCAID5yHwCAAzQMeAIAABmJfgEAGAxmwafVXVDVT1SVfeutP1yVf15Vd1TVR+oqmdN7edV1V9X1d3T8vaV57yoqj5RVbur6q1VtX55ZgCAAzpJ17hlkKPJfN6Y5LKD2m5L8l3d/d1J/mOSN65s+3R3XzgtP7nSfl2SH09y/rQcfEwAAHa4TYPP7v5okkcPavuD7n5ienh7knOOdIyqOivJM7v79u7uJDclecXxdRkAYGfoHreMciK+8/kvk/zeyuPnVtWfVdW/r6rvn9rOTrJnZZ89UxsAAGtkS7daqqp/k+SJJO+amh5K8pzu/mJVvSjJB6vq+cdx3F1JdiXJ0/INW+kiAADbyHEHn1X16iQ/nOTFUyk93f14ksen9Tur6tNJnpdkb76+NH/O1HZI3X19kuuT5Jn17DW89z8AsBbWMMo5rrJ7VV2W5F8neXl3f2Wl/Zur6pRp/duzMbDoM939UJLHquqSaZT7q5J8aMu9BwDgpLJp5rOqbk5yaZIzqmpPkmuzMbr9tCS3TXdMun0a2f4DSd5UVV9Nsj/JT3b3gcFKP5WNkfN/NxvfEV39nigAwJqptZzhaNPgs7uvOkTzOw6z7/uSvO8w2+5I8l3H1DsAAHYUMxwBALCYLY12BwBgCww4AgCA+Qg+AQBYjLI7AMAInbUc7S7zCQDAYmQ+AQBGMeAIAADmI/gEAGAxyu4AAMMYcAQAALOR+QQAGMWAIwAAmI/gEwCAxQg+AQBG6YHLJqrqhqp6pKruPcz2S6vqS1V197T83NGcsu98AgBwKDcmeVuSm46wzx939w8fy0EFnwAAI3SSbTy3e3d/tKrOO9HHVXYHAOB4/cOq+g9V9XtV9fyjeYLMJwDAejqjqu5YeXx9d19/DM+/K8m3dfdfVdXlST6Y5PzNniT4BAAYpMfe5/ML3X3R8T65ux9bWb+1qn6jqs7o7i8c6XnK7gAAHLOq+paqqmn94mzElV/c7HkynwAAo2zjGY6q6uYkl2ajPL8nybVJTk2S7n57kh9N8q+q6okkf53kyu7Nc7mCTwAAnqS7r9pk+9uycSumY6LsDgDAYmQ+AQBG2cb3+ZyLzCcAAIsRfAIAsBhldwCAQWobj3afi8wnAACLkfkEABihs63v8zkXmU8AABYj+AQAYDHK7gAAQ5T7fAIAwJxkPgEARjHgCAAA5iP4BABgMcruAACjKLsDAMB8ZD4BAEaR+QQAgPkIPgEAWIyyOwDACB0zHAEAwJxkPgEABikDjgAAYD6CTwAAFqPsDgAwirI7AADMR/AJAMBiBJ8AACxG8AkAwGIMOAIAGMR9PgEAYEaCTwAAFqPsDgAwStfoHixO5hMAgMUcVfBZVTdU1SNVde9K27Or6raq+tT07+lTe1XVW6tqd1XdU1UvXHnO1dP+n6qqq0/86QAAnCR68DLI0WY+b0xy2UFtb0jyke4+P8lHpsdJ8rIk50/LriTXJRvBapJrk3xPkouTXHsgYAUAYD0cVfDZ3R9N8uhBzVckeee0/s4kr1hpv6k33J7kWVV1VpKXJrmtux/t7r9MclueHNACALCDbWXA0Znd/dC0/vkkZ07rZyd5cGW/PVPb4dqfpKp2ZSNrmqflG7bQRQCAbcx9Po9Pd5/Qbw909/XdfVF3X3RqTjtRhwUAYLCtBJ8PT+X0TP8+MrXvTXLuyn7nTG2HawcAWEvV45ZRthJ83pLkwIj1q5N8aKX9VdOo90uSfGkqz384yUuq6vRpoNFLpjYAANbEUX3ns6puTnJpkjOqak82Rq3/YpL3VtVrk3wuyY9Nu9+a5PIku5N8JclrkqS7H62qX0jy8Wm/N3X3wYOYAADYwY4q+Ozuqw6z6cWH2LeTXHOY49yQ5Iaj7h0AwE5mwBEAAMzH3O4AAKPIfAIAwHwEnwAALEbZHQBggNH32xxF5hMAgMXIfAIAjNI1ugeLk/kEAGAxgk8AABaj7A4AMIoBRwAAMB+ZTwCAQdxqCQAAZiT4BABgMcruAACjKLsDAMB8BJ8AACxG2R0AYIQ22h0AAGYl8wkAMIrMJwAAzEfwCQDAYpTdAQBGUXYHAID5yHwCAAziVksAADAjwScAAIsRfAIAsBjBJwAAizHgCABgFAOOAABgPoJPAAAWo+wOADBCu88nAADMSuYTAGAUmU8AAJiP4BMAgMUouwMAjKLsDgAA8xF8AgCwGGV3AIABKu7zCQAAs5L5BAAYReYTAADmI/gEAGAxyu4AACO0AUcAADArmU8AgFFkPgEAYD6CTwAAFqPsDgAwirI7AADMR+YTAGAQt1oCAIAZCT4BAFiMsjsAwCjK7kevqr6zqu5eWR6rqtdX1c9X1d6V9stXnvPGqtpdVZ+sqpeemFMAAOBkcdyZz+7+ZJILk6SqTkmyN8kHkrwmya9196+s7l9VFyS5Msnzk3xrkj+squd1977j7QMAwEmrI/O5BS9O8unu/twR9rkiybu7+/Hu/osku5NcfIJeHwCAk8CJCj6vTHLzyuPXVdU9VXVDVZ0+tZ2d5MGVffZMbQAArIktB59V9dQkL0/yO1PTdUm+Ixsl+YeSvOU4jrmrqu6oqju+mse32kUAgG2petwyyonIfL4syV3d/XCSdPfD3b2vu/cn+c18rbS+N8m5K887Z2p7ku6+vrsv6u6LTs1pJ6CLAABsByci+LwqKyX3qjprZduPJLl3Wr8lyZVVdVpVPTfJ+Un+9AS8PgDAyakHLoNs6T6fVfX0JD+U5CdWmn+pqi7Mxml99sC27r6vqt6b5P4kTyS5xkh3AID1sqXgs7v/c5K/d1DbK4+w/5uTvHkrrwkAwMnLDEcAAIOMHPgzirndAQBYjOATAIDFKLsDAIyi7A4AAMk0U+UjVXXvYbZXVb21qnZPM1u+8GiOK/gEABhh5D0+jy7jemOSy46w/WXZuG/7+Ul2ZWOWy00JPgEAeJLu/miSR4+wyxVJbuoNtyd51kGTDR2S4BMAYD2dUVV3rCy7jvH5Zyd5cOXxnqntiAw4AgAYoKZloC9090VLv6jMJwAAx2NvknNXHp8ztR2R4BMAYJTtPeBoM7ckedU06v2SJF/q7oc2e5KyOwAAT1JVNye5NBvfDd2T5NokpyZJd789ya1JLk+yO8lXkrzmaI4r+AQA4Em6+6pNtneSa471uIJPAIBBygxHAAAwH5lPAIBRZD4BAGA+gk8AABaj7A4AMIqyOwAAzEfmEwBghHarJQAAmJXgEwCAxSi7AwCMouwOAADzkfkEABjEgCMAAJiR4BMAgMUouwMAjKLsDgAA8xF8AgCwGGV3AIBBjHYHAIAZyXwCAIzQMeAIAADmJPgEAGAxyu4AAKMouwMAwHxkPgEABqi41RIAAMxK8AkAwGKU3QEARlF2BwCA+ch8AgAMUr1+qU+ZTwAAFiP4BABgMcruAAAjdAw4AgCAOcl8AgAMYoYjAACYkeATAIDFKLsDAIyi7A4AAPMRfAIAsBhldwCAQYx2BwCAGW05+Kyqz1bVJ6rq7qq6Y2p7dlXdVlWfmv49fWqvqnprVe2uqnuq6oVbfX0AgJNWD1wGOVGZzx/s7gu7+6Lp8RuSfKS7z0/ykelxkrwsyfnTsivJdSfo9QEAOAnMVXa/Isk7p/V3JnnFSvtNveH2JM+qqrNm6gMAANvMiQg+O8kfVNWdVbVrajuzux+a1j+f5Mxp/ewkD648d8/UBgCwXnpjwNGoZZQTMdr9H3X33qr6+0luq6o/X93Y3V11bKc4BbG7kuRp+YYT0EUAALaDLWc+u3vv9O8jST6Q5OIkDx8op0//PjLtvjfJuStPP2dqO/iY13f3Rd190ak5batdBADYngw4OjZV9fSq+sYD60lekuTeJLckuXra7eokH5rWb0nyqmnU+yVJvrRSngcAYIfbatn9zCQfqKoDx/rt7v79qvp4kvdW1WuTfC7Jj03735rk8iS7k3wlyWu2+PoAAJxEthR8dvdnkvyDQ7R/McmLD9HeSa7ZymsCAOwEFTMcAQDArMztDgAwSq9f6lPmEwCAxQg+AQBYjLI7AMAgBhwBAMCMZD4BAEYYPNPQKDKfAAAsRvAJAMBilN0BAAap/aN7sDyZTwAAFiPzCQAwigFHAAAwH8EnAACLUXYHABjEDEcAADAjwScAAItRdgcAGKGT9PrV3WU+AQBYjMwnAMAgBhwBAMCMBJ8AACxG2R0AYBRldwAAmI/MJwDAABUDjgAAYFaCTwAAFqPsDgAwQrcZjgAAYE4ynwAAgxhwBAAAMxJ8AgCwGGV3AIBRlN0BAGA+Mp8AAIMYcAQAADMSfAIAsBhldwCAETrJ/vWru8t8AgCwGJlPAIBR1i/xKfMJAMByBJ8AACxG2R0AYBD3+QQAgBkJPgEAWIyyOwDAKL1+dXeZTwAAFiPzCQAwiAFHAAAwI8EnAACLUXYHABihY3pNAACYk8wnAMAAlaTcagkAAOYj+AQAYDHK7gAAo+wf3YHlyXwCALCY4w4+q+rcqvqjqrq/qu6rqp+e2n++qvZW1d3TcvnKc95YVbur6pNV9dITcQIAACer6h62jLKVsvsTSX6mu++qqm9McmdV3TZt+7Xu/pXVnavqgiRXJnl+km9N8odV9bzu3reFPgAAcBI57sxndz/U3XdN619O8kCSs4/wlCuSvLu7H+/uv0iyO8nFx/v6AADMq6oumyrWu6vqDYfY/uqq+v9WKt7/w2bHPCHf+ayq85K8IMnHpqbXVdU9VXVDVZ0+tZ2d5MGVp+3JYYLVqtpVVXdU1R1fzeMnoosAANtLD142UVWnJPn1JC9LckGSq6ZK9sHe090XTstvbXbcLQefVfWMJO9L8vrufizJdUm+I8mFSR5K8pZjPWZ3X9/dF3X3RafmtK12EQCAY3dxkt3d/Znu/psk785GJXtLthR8VtWp2Qg839Xd70+S7n64u/d19/4kv5mvldb3Jjl35ennTG0AAGuokx64bO5oq9b/bKp4/25VnXuI7V9nK6PdK8k7kjzQ3b+60n7Wym4/kuTeaf2WJFdW1WlV9dwk5yf50+N9fQAAtuSMA19znJZdx3GMf5vkvO7+7iS3JXnnZk/Yymj370vyyiSfqKq7p7afzcb3AS7MxrcJPpvkJ5Kku++rqvcmuT8bI+WvMdIdAGCYL3T3RUfYvmnVuru/uPLwt5L80mYvetzBZ3f/SZI6xKZbj/CcNyd58/G+JgDATlLjbrd5ND6e5PypYr03G7fM/O9Wd6iqs7r7oenhy7Nx96MjMr0mAABP0t1PVNXrknw4ySlJbpgq2W9Kckd335Lkf6qql2ejqv1okldvdlzBJwAAh9Tdt+agqnZ3/9zK+huTvPFYjin4BAAYZeA0l6OckJvMAwDA0ZD5BAAYoZPaP7oTy5P5BABgMYJPAAAWo+wOADCKAUcAADAfmU8AgFHWL/Ep8wkAwHIEnwAALEbZHQBgkDLgCAAA5iPzCQAwiswnAADMR/AJAMBilN0BAEboJPtHd2J5Mp8AACxG5hMAYIBKu9USAADMSfAJAMBilN0BAEZRdgcAgPnIfAIAjKrdAykAAAnLSURBVCLzCQAA8xF8AgCwGGV3AIARzHAEAADzEnwCALAYZXcAgEFMrwkAADOS+QQAGEXmEwAA5iP4BABgMcruAABDtLI7AADMSeYTAGCEjswnAADMSfAJAMBilN0BAEbZP7oDy5P5BABgMTKfAACDmNsdAABmJPgEAGAxyu4AAKMouwMAwHxkPgEARugk+2U+AQBgNoJPAAAWo+wOADBEG3AEAABzkvkEABhF5hMAAOYj+AQAYDHK7gAAoyi7z6+qLquqT1bV7qp6w9KvDwDAOIsGn1V1SpJfT/KyJBckuaqqLliyDwAAjLN02f3iJLu7+zNJUlXvTnJFkvsX7gcAwFim11zE2UkeXHm8Z2oDAGANbMsBR1W1K8mu6eHjf9i/e+/I/mwDZyT5wuhODOYabHAdXIMDXAfXIHENDjjW6/Btc3Xk2HTS+0d3YnFLB597k5y78vicqe3rdPf1Sa5Pkqq6o7svWqZ725Nr4Boc4Dq4Bge4Dq5B4hoc4DqcXJYuu388yflV9dyqemqSK5PcsnAfAAAYZNHMZ3c/UVWvS/LhJKckuaG771uyDwAA28Ya3udz8e98dvetSW49hqdcP1dfTiKugWtwgOvgGhzgOrgGiWtwgOtwEqlew4gbAGC0b3rqmf2933LVsNf//Qf/9ztHfFfW3O4AACxm2waf6zINZ1WdW1V/VFX3V9V9VfXTU/vPV9Xeqrp7Wi5fec4bp+vyyap66bjen1hV9dmq+sR0vndMbc+uqtuq6lPTv6dP7VVVb52uwz1V9cKxvd+6qvrOlZ/33VX1WFW9fh3eC1V1Q1U9UlX3rrQd88++qq6e9v9UVV094lyO12GuwS9X1Z9P5/mBqnrW1H5eVf31ynvi7SvPedH0e7R7uk414nyO12GuwzH/DpzMnyGHuQbvWTn/z1bV3VP7jnwvHOGzca3+LuxU2zL4rPWahvOJJD/T3RckuSTJNSvn+mvdfeG03Jok07Yrkzw/yWVJfmO6XjvFD07ne6AM8IYkH+nu85N8ZHqcbLw3zp+WXUmuW7ynJ1h3f/LAzzvJi5J8JckHps07/b1wYzbOYdUx/eyr6tlJrk3yPdmYTe3aAx9MJ4kb8+RrcFuS7+ru707yH5O8cWXbp1feEz+50n5dkh/P167Rwcfc7m7Moft81L8DO+Az5MYcdA26+1+s/H14X5L3r2zeie+Fw3027ry/C93jlkG2ZfCZlWk4u/tvkhyYhnPH6e6Huvuuaf3LSR7IkWd9uiLJu7v78e7+iyS7s3G9dqorkrxzWn9nklestN/UG25P8qyqOmtEB2fy4mx8oHzuCPvsmPdCd380yaMHNR/rz/6lSW7r7ke7+y+zEbidNB+2h7oG3f0H3f3E9PD2bNwb+bCm6/DM7r69N77Qf1O+dt1OCod5LxzO4X4HTurPkCNdgyl7+WNJbj7SMU7298IRPhvX6u/CTrVdg8+1nIazqs5L8oIkH5uaXjeVD25Y+Z/aTr42neQPqurO2pjlKknO7O6HpvXPJzlzWt/J1yHZyOasfris23shOfaf/U6/Hv8yye+tPH5uVf1ZVf37qvr+qe3sbJz3ATvpGhzL78BOfi98f5KHu/tTK207+r1w0Gfjzvu7IPPJKFX1jGyUUl7f3Y9lo2TwHUkuTPJQkrcM7N5S/lF3vzAb5ZNrquoHVjdO/3vf8bdnqI0JGF6e5HempnV8L3yddfnZH05V/ZtslCHfNTU9lOQ53f2CJP9zkt+uqmeO6t8C1v53YMVV+fr/mO7o98IhPhv/1rr/XTiZbdfg86im4dwpqurUbPxyvau7358k3f1wd+/r7v1JfjNfK6fu2GvT3Xunfx/JxncdL07y8IFy+vTvI9PuO/Y6ZCP4vqu7H07W870wOdaf/Y68HlX16iQ/nOS/nz5sM5WZvzit35nk00mel43zXS3N74hrcBy/Azv1vfCUJP80yXsOtO3k98KhPhvj78KOsF2Dz7WZhnP6/s47kjzQ3b+60r76/cUfSXJg1OMtSa6sqtOq6rnZ+HL1ny7V37lU1dOr6hsPrCd5STbO+ZYkB0YnXp3kQ9P6LUleNY1wvCTJl1ZKMSe7r8tsrNt7YcWx/uw/nOQlVXX6VJZ9ydR20qqqy5L86yQv7+6vrLR/84HBZVX17dn42X9mug6PVdUl09+WV+Vr1+2kdRy/Azv1M+SfJPnz7v7bcvpOfS8c7rMxO+7vwsCS+8Cy++IzHB2NNZuG8/uSvDLJJ2q6dUaSn83G6MwLs1FS+GySn0iS7r6vqt6b5P5slOGu6e59i/f6xDszyQc2/t7kKUl+u7t/v6o+nuS9VfXaJJ/Lxhftk41Zsi7PxgCDryR5zfJdPvGmwPuHMv28J7+0098LVXVzkkuTnFFVe7IxOvUXcww/++5+tKp+IRuBR5K8qbuPduDKcIe5Bm9MclqS26bfjdun0cw/kORNVfXVJPuT/OTKuf5UNkZL/91sfEd09Xui295hrsOlx/o7cDJ/hhzqGnT3O/Lk74InO/e9cLjPxrX6u7BTmeEIAGCAbzr17/f3nvHPh73+73/+N8xwBADAzib4BABgMdvyO58AAGthDb/+KPMJAMBiZD4BAEaR+QQAgPkIPgEAWIyyOwDAEJ3sV3YHAIDZCD4BAFiMsjsAwAiddO8f3YvFyXwCALAYmU8AgFEMOAIAgPkIPgEAWIyyOwDAKKbXBACA+ch8AgCM0J3sd6slAACYjeATAIDFKLsDAIxiwBEAAMxH5hMAYJA24AgAAOYj+AQAYDHK7gAAQ7QBRwAAMCeZTwCAETrJfplPAACYjeATAIDFKLsDAIzS7vMJAACzEXwCALAYZXcAgAE6SRvtDgAA85H5BAAYoduAIwAAmJPgEwCAxQg+AQAG6f09bDkaVXVZVX2yqnZX1RsOsf20qnrPtP1jVXXeZscUfAIA8CRVdUqSX0/ysiQXJLmqqi44aLfXJvnL7v6vkvxakv9ts+MKPgEARun945bNXZxkd3d/prv/Jsm7k1xx0D5XJHnntP67SV5cVXWkgwo+AQA4lLOTPLjyeM/Udsh9uvuJJF9K8veOdFC3WgIAGODL+csP/2H/7hkDu/C0qrpj5fH13X393C8q+AQAGKC7Lxvdh03sTXLuyuNzprZD7bOnqp6S5JuSfPFIB1V2BwDgUD6e5Pyqem5VPTXJlUluOWifW5JcPa3/aJJ/191HHEov8wkAwJN09xNV9bokH05ySpIbuvu+qnpTkju6+5Yk70jyf1bV7iSPZiNAPaLaJDgFAIATRtkdAIDFCD4BAFiM4BMAgMUIPgEAWIzgEwCAxQg+AQBYjOATAIDFCD4BAFjM/w+bjQszQFxdtwAAAABJRU5ErkJggg==\n",
      "text/plain": [
       "<Figure size 864x864 with 2 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plot_image(output_model.data, index=(0,0))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.7.4"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
