OLDVERSION=v1.08
VERSION=v1.09

TARGET=/Users/gaser/spm/spm2/toolbox/vbm2

STARGET=141.35.200.101:/Applications/xampp/htdocs/

FILES=cg_create_template.m cg_gwc_HMRF.m cg_read_vbm_volumes.m cg_spm_ui.m cg_vbm_longitudinal_bias.m cg_vbm_optimized.m cg_vbm_lat_index.m cg_check_sample_sd.m cg_spmT2x.m cg_showslice_all.m spm_segment.m spm_segment_ui.m spm_list.m spm_getSPM.m spm_normalise.m spm_max_nS.m spm_spm.m spm_vbm2.m vbm2.man INSTALL.txt

ZIPFILE=vbm2_$(VERSION).zip

install:
	-@test ! -d $(TARGET) || rm -r $(TARGET)
	-@mkdir $(TARGET)
	-@cp $(FILES) $(TARGET)

upgrade:
	-@for i in $(FILES); do sed -i "" -e "s/$(OLDVERSION)/$(VERSION)/g" $$i; done    

help:
	-@echo Available commands:
	-@echo install zip scp upgrade

zip: upgrade
	-@test ! -d vbm2 || rm -r vbm2
	-@cp -r $(TARGET) .
	-@zip $(ZIPFILE) -rm vbm2

scp:    zip
	-@scp -pr $(ZIPFILE) $(STARGET)
