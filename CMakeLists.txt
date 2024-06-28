add_library(fabm_models_niva OBJECT
            niva_model_library.F90
            brom/brom_acetate.F90
            brom/brom_bact.F90
            brom/brom_bio.F90
            brom/brom_calcium.F90
            brom/brom_carbon.F90
            brom/brom_cod.F90
            brom/brom_eq_constants.F90
            brom/brom_fe.F90
            brom/brom_hg.F90
            brom/brom_main_nutrients.F90
            brom/brom_manganese.F90
            brom/brom_methane.F90
            brom/brom_ni.F90
            brom/brom_nitrogen.F90
            brom/brom_partitioning.F90
            brom/brom_pH.F90
            brom/brom_salt.F90
            brom/brom_silicon.F90
            brom/brom_sulfur.F90
            brom/brom_ba.F90
            brom/brom_bubble.F90
            brom/oxydep.F90
            brom/brom_injection.F90
            brom/brom_bioplast.F90
	        brom/brom_mg.F90
            brom/light.F90
           )

add_dependencies(fabm_models_niva fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

