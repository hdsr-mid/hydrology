Updates:


07-08-2024:
- functies_droog_profie.py	coordinaten afgerond naar 2 decimalen in variabele ind_drop (vanaf lijn 346)
- functies_droog_profie.py	d_edge opnieuw berekend voor iedere conditie in lijn 427 (anders wordt het namelijk overschreven van de vorige iteratie)
- functies_droog_profie.py	variabele d_extra verwijderd zodat extra punten op bepaalde afstanden berekend worden

23-09-2024
- functies_droog_profie.py	nieuwe variabele d_extrapolate = 15 voor het extrapoleren van profielpunten voorbij de insteeklijn. Eerder werd 5 m als grens gehanteerd welke te klein bleek te zijn. Deze variabele is gebruikt in lijnen 148 - 152
- functies_droog_profie.py	continue -> break (lijn 560 en 591); zodat er alleen extra punten toegevoegd worden zolang de hoogte toeneemt. Zodra de hoogte afneemt, worden geen punten meer toegevoegd.