Doel: DHydro model clippen op afvoergebied

Basismodel: Hertoetsing, opgeleverd sep 2025

Stappen:
- Clip branches-shapefile o.b.v. afvoergebied-polygoon (shapefile)
- Clip dflowfm/network.nc o.b.v. geclipte branches
- Clip bijbehorende kuntwerken, profiel info, laterals, observation points, boundary conditions, initial conditions
- Clip bedlevel.xyz en hoogtelijnen (.pliz) o.b.v. box-coordinaten afgeleid van het afvoergebied-polygoon (shapefile)
- Clip tif o.b.v. polygon
- Clip RR-knopen
- Clip RTC-knopen/sturing
- Verwijder overbodige koppelingen in dimr
- HANDMATIG in GUI: 
    - pas boundary conditions aan
    - controleer RTC-sturing (ontbreekt er iets?)


Aanname:
- er is nu een buffer van 1m rondom het grensgebied gehanteerd


Een paar codes...
clip_code      = 'AFVGEB0057'
naam           = 'Molenvliet' 
    
clip_code      = 'AFVGEB0035'
naam           = 'KrommeRijn'

clip_code      = 'AFVGEB0040'
naam           = 'ARK'

clip_code      = 'AFVGEB0055'
naam           = 'EvS'

clip_code      = 'AFVGEB0100' # IJssel
naam           = 'IJssel'

clip_code      = 'AFVGEB0052'
naam           = 'OudeRijn'

clip_code      = 'AFVGEB0097'
naam           = 'LeidseRijn'

clip_code      = 'AFVGEB0092'
naam           = 'Honswijk'

clip_code      = 'AFVGEB0064'
naam           = 'Heemstede'

clip_code      = 'AFVGEB0065'
naam           = 'Poeldijk'