Contact:
Petra Husman

Hoofdscript van waaruit gerund kan worden:
1_Workflow_model_verificatie.ipynb
2_Workflow_model_validatie.ipynb
3_Workflow_model_evaluatie.ipynb

Benodigde virtuele omgeving:
Petra_arcpy (op reken05)
env_merge_html (op vdi)
geo-env + extras (MaaS server)

Doel van iedere subscript:
- functies_punt_laterals_T1.py		3	validatie afvoer (laterals)			RUPROF input vs. Hydromedah (via Sobek model)
- functies_punt_Q_T1.py			2	validatie debiet, T=1, zomer			RUPROF vs. FEWS-WIS
- functies_punt_Q_T1_bokeh.py		2	validatie debiet, T=1, zomer			
- functies_punt_Qdir.py			2	validatie stroomrichting			Sobek model vs. satelliet data
- functies_punt_Qdir_bokeh.py		2	validatie stroomrichting, plotting		
- functies_punt_Qseries.py		2	validatie debiet				Sobek model vs. satelliet data
- functies_punt_Qseries_bokeh.py	2	validatie debiet, plotting
- functies_punt_WL_T1.py		2	validatie waterniveau, T=1, zomer		RUPROF vs. FEWS-WIS
- functies_punt_WL_T1_bokeh.py		2	validatie waterniveau, T=1, zomer
- functies_ruimtelijk_afvoer_gem.py	3	validatie afvoer (laterals)			Sobek model vs. satelliet data
- functies_ruimtelijk_peil.py		1	validatie peil					BR streefpeil vs. peilbesluit
- functies_ruimtelijk_streefpeil.py	1	validatie streefpeil				RUPROF waterniveau op t=0 vs. BR streefpeil 