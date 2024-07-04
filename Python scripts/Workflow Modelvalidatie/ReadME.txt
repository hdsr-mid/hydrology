Contact:
Petra Husman

Hoofdscript van waaruit alles gerund kan worden:
Workflow_model_validatie.ipynb

Benodigde virtuele omgeving:
Petra_arcpy

Doel van iedere subscript:
- functies_punt_laterals_T1.py			validatie afvoer (laterals)			RUPROF input vs. Hydromedah (via Sobek model)
- functies_punt_Q_T1.py				validatie debiet, T=1, zomer			RUPROF vs. FEWS-WIS
- functies_punt_Qdir.py				validatie stroomrichting			Sobek model vs. satelliet data
- functies_punt_Qdir_bokeh.py			validatie stroomrichting, plotting		
- functies_punt_Qseries.py			validatie debiet				Sobek model vs. satelliet data
- functies_punt_Qseries_bokeh.py		validatie debiet, plotting
- functies_punt_WL_T1.py			validatie waterniveau, T=1, zomer		RUPROF vs. FEWS-WIS
- functies_ruimtelijk_afvoer_gem.py		validatie afvoer (laterals)			Sobek model vs. satelliet data
- functies_ruimtelijk_streefpeil.py		validatie streefpeil				RUPROF waterniveau op t=0 vs. BR streefpeil 