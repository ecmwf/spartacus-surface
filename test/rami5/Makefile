# This Makefile is for running SPARTACUS-Surface on the RAMI-V
# "actual" test scenes.  The variables required by SPARTACUS are
# provided in the NetCDF files in the scene_nc directory.  Typing
# "make" will run SPARTACUS on them to produce outputs in the out_nc
# directory. The Matlab script "process_spartacus_scenes.m" produces
# the ASCII files for submission in the "mes" directory.

BASE = rami5
SPSURF = ../../bin/spartacus_surface
CONFIG = config.nam
CONFTMP = config_tmp.nam
CHANGE_NL = ./change_namelist.sh
INDIR = scene_nc
OUTDIR = out_nc

INSUFFIX = _scene.nc
OUTSUFFIX = _out.nc
LOGSUFFIX = _out.log
BSOUTSUFFIX = _blacksoil_out.nc
BSLOGSUFFIX = _blacksoil_out.log

# Default target when you run "make": run the five actual actual cases
# with solar configurations required by RAMI-V
test: rami5
rami5: HET09_JBS_SUM HET07_JPS_SUM HET08_OPS_WIN HET14_WCO_UND HET15_JBS_WIN

# For a larger number of runs when you run "make all": also includes
# overhead sun and the RAMI-VI sun angles
all: rami5 overhead rami4

# Individual targets consisting of <SCENE>-<ILLUMINATION> where
# <ILLUMINATION> is either "diffuse" indicating the incoming solar
# radiation at top-of-canopy is isotropic, or <SZA>-direct where SZA
# is the solar zenith angle in degrees

# Jarvselja Birch Stand - Summer
HET09_JBS_SUM: HET09_JBS_SUM-diffuse HET09_JBS_SUM-56-direct HET09_JBS_SUM-41-direct
# Jarvselja Birch Stand - Winter
HET15_JBS_WIN: HET15_JBS_WIN-diffuse HET15_JBS_WIN-76-direct HET15_JBS_WIN-56-direct
# Jarvselja Pine Stand - Summer
HET07_JPS_SUM: HET07_JPS_SUM-diffuse HET07_JPS_SUM-56-direct HET07_JPS_SUM-41-direct
# Offenpass Pine Stand - Winter
HET08_OPS_WIN: HET08_OPS_WIN-diffuse HET08_OPS_WIN-76-direct HET08_OPS_WIN-56-direct
# Wellington Citrus Orchard
HET14_WCO_UND: HET14_WCO_UND-diffuse HET14_WCO_UND-42-direct HET14_WCO_UND-60-direct HET14_WCO_UND-67-direct

# Run the five scenes with overhead sun, from which the penetration of
# direct sunlight is computed and can be compared with values on the
# RAMI-V web site
overhead: HET09_JBS_SUM-00-direct HET15_JBS_WIN-00-direct  HET07_JPS_SUM-00-direct HET08_OPS_WIN-00-direct HET14_WCO_UND-00-direct

# Solar zenith angles use by the same scenes from the RAMI-IV
# intercomparison
rami4: HET09_JBS_SUM-37-direct HET15_JBS_WIN-54-direct HET07_JPS_SUM-37-direct HET08_OPS_WIN-47-direct HET14_WCO_UND-00-direct HET14_WCO_UND-20-direct HET14_WCO_UND-50-direct

# Run SPARTACUS-Surface on scene with diffuse top-of-canopy solar
# illumination, both with the real soil albedo and a "black soil"
# configuration
%-diffuse:
	mkdir -p out_nc
	set -o pipefail && $(SPSURF) $(CONFIG) $(INDIR)/$(BASE)_$*$(INSUFFIX) \
		$(OUTDIR)/$(BASE)_$*-diffuse$(OUTSUFFIX) | tee $(OUTDIR)/$(BASE)_$*-diffuse$(LOGSUFFIX)
	$(CHANGE_NL) $(CONFIG) $(CONFTMP) ground_sw_albedo=0.0
	set -o pipefail && $(SPSURF) $(CONFTMP) $(INDIR)/$(BASE)_$*$(INSUFFIX) \
		$(OUTDIR)/$(BASE)_$*-diffuse$(BSOUTSUFFIX) | tee $(OUTDIR)/$(BASE)_$*-diffuse$(BSLOGSUFFIX)

# Run SPARTACUS-Surface on scene with direct top-of-canopy solar
# illumination.  First create config_tmp.nam with the required solar
# zenith angle, then run SPARTACUS-Surface, in both cases using awk to
# extract parts of the "stem" represented by the "%" character in GNU
# make.
%-direct:
	mkdir -p out_nc
	$(CHANGE_NL) $(CONFIG) $(CONFTMP) top_flux_dn_direct_sw=1.0 \
		solar_zenith_angle=$(shell echo $* | awk -F- '{print $$2}')
	set -o pipefail && $(SPSURF) $(CONFTMP) $(INDIR)/$(BASE)_$(shell echo $* | awk -F- '{print $$1}')$(INSUFFIX) \
		$(OUTDIR)/$(BASE)_$*$(OUTSUFFIX) | tee $(OUTDIR)/$(BASE)_$*$(LOGSUFFIX)
	$(CHANGE_NL) $(CONFIG) $(CONFTMP) top_flux_dn_direct_sw=1.0 \
		solar_zenith_angle=$(shell echo $* | awk -F- '{print $$2}') \
		ground_sw_albedo=0.0
	set -o pipefail && $(SPSURF) $(CONFTMP) $(INDIR)/$(BASE)_$(shell echo $* | awk -F- '{print $$1}')$(INSUFFIX) \
		$(OUTDIR)/$(BASE)_$*$(BSOUTSUFFIX) | tee $(OUTDIR)/$(BASE)_$*$(BSLOGSUFFIX)

clean:
	rm -rf $(CONFTMP) out_nc/$(BASE)*.nc out_nc/$(BASE)*.log mes/*.mes mes/*.png
