Workflow:
- Clip AHN raster
- Bewerk AHN raster (gaten opvullen, hoogte watervlaktes gelijk zetten aan streefpijl)
- Sobek waterstand data aan shapefile koppelen (nodes)
- Max waterstand per peilgebied bepalen
- Inundatie berekenen o.b.v. max waterstand voor ieder peilgebied
- Inundatiekaarten per peilgebied met elkaar mergen
- Inundatiegebieden "verwijderen" ("als droog markeren") als deze gridcellen niet (via andere natte gridcellen) verbonden zijn met een watergang 
(-> een grid cell kan lager liggen dan de max waterstand, maar omringd door hoger gelegen gebieden waardoor deze alsnog droog blijft tenzij er een duiker is oid)
- Inundatiekaarten plotten