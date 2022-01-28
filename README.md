# Opérateurs d'intérêt

### Détection des contours
Toute l’idée de détection des points isolés et coins repose sur la détection des contours. Là où les contours sont à plusieurs directions, on considère le point comme étant un coin : Détecter pixels où le gradient est fort dans plus d’une direction – ce sont les pixels situés sur les coins.
