Updates:

30-01-2025:
- functies_droog_profiel.py	Datum WIT meting meenemen in de output shapefile (lijnen 28-29)

16-12-2024:
- functies_droog_profiel.py	Nieuwe filter: Profielen worden verwijderd als hun lijnen uit meerdere stukken bestaan (teken dat het meerdere insteekvlakken dekt). Lijnen 267 - 275 zijn toegevoegd.

14-11-2024:
- functies_droog_profiel.py	Nieuwe filter: De hydro object code van het meetpunt het dichtst bij een hydro object moet gelijk zijn aan de code van het insteekvlak. Im sommige gevallen ligt het profiel namelijk net bij een ander hydroobject (voorbeeld profiel 1896076). Dit is toegevoegd in lijnen 56 - 65 en 98 - 106
- functies_droog_profiel.py	De lengte van de ge-extrapoleerde lijn was oorspronkelijk ook afhankelijk van de breedte van het gemeten profiel. In sommige gevallen was deze breedte significant kleiner dan de breedte o.b.v. de insteek. Daarom is dit verwijderd in lijnen 173 - 177 (voorbeeld profiel 1893856).
- functies_droog_profiel.py	Lijn 266 is toegevoegd zodat een ge-extrapoleerde lijn alleen 1 watergang kruist. Soms liggen namelijk twee watergangen erg dicht bij elkaar (voorbeeld profiel 1899641).
- functies_droog_profiel.py	Lijnen 270 - 272 zijn verwijderd omdat het niet nodig is die variabelen op te slaan (gdf_single, gdf_multi, gdf_multi_joined).

23-09-2024
- functies_droog_profiel.py	nieuwe variabele d_extrapolate = 15 voor het extrapoleren van profielpunten voorbij de insteeklijn. Eerder werd 5 m als grens gehanteerd welke te klein bleek te zijn. Deze variabele is gebruikt in lijnen 148 - 152
- functies_droog_profiel.py	continue -> break (lijn 560 en 591); zodat er alleen extra punten toegevoegd worden zolang de hoogte toeneemt. Zodra de hoogte afneemt, worden geen punten meer toegevoegd.

07-08-2024:
- functies_droog_profiel.py	coordinaten afgerond naar 2 decimalen in variabele ind_drop (vanaf lijn 346)
- functies_droog_profiel.py	d_edge opnieuw berekend voor iedere conditie in lijn 427 (anders wordt het namelijk overschreven van de vorige iteratie)
- functies_droog_profiel.py	variabele d_extra verwijderd zodat extra punten op bepaalde afstanden berekend worden

