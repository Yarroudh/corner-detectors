# Opérateurs d'intérêt

Détection des contours
Toute l’idée de détection des points isolés et coins repose sur la détection des contours. Là où les contours sont à plusieurs directions, on considère le point comme étant un coin.
	Détecter pixels où le gradient est fort dans plus d’une direction – ce sont les pixels situés sur les coins
 
Méthode de Harris
Etape 1 : Calculer les images de gradient I_x et I_y 
	La détection des contours peut se faire facilement à l’aide de l’opérateur Sobel. Il y a d’autres opérateurs 	comme celui de Prewitt, celui que nous avant utilisé pour calculer le gradient de la pente à partir d’un MNT.
	Opérateur de Sobel :
I_x=[■(-1&0&1@-2&0&2@-1&0&1)]×A 		I_y=[■(-1&-2&-1@0&0&0@1&2&1)]×A

Etape 2 : pour tous les points de l’image de gradient, calculer la matrice de covariance du gradient (dans une fenêtre 2N+1 × 2N+1) et en extraire les valeurs propres λ_1 et λ_2
	Pour visualiser cette étape, on peut prendre un pixel (x,y) et un voisinage (ΔX,ΔY)
 
	La fonction d’autocorrélation de l’image en (x,y) est donnée par :
	Intuitivement, la corrélation de deux objets mesure leur dépendance réciproque ; L'autocorrélation d'un 	signal mesurera donc les dépendances internes de ce signal. Nous l'appliquerons dans le cas d'un signal 	à une variable entière.
	Une application sur laquelle nous nous attarderons un peu ici est la détection du bruit dans une 	image. 
c(x,y)=∑_W▒〖[I(x_i,y_i )-I(〗 x_i+Δx,y_i+Δy)]²
	Si on approxime I(x_i+Δx,y_i+Δy) par sa série de Taylor, en se limitant aux 2 premiers termes, on 	obtient :
I(x_i+Δx,y_i+Δy)=I(x_i,y_i )+∂I/∂x |_(x_i,y_i ) Δx+∂I/∂y |_(x_i,y_i ) Δy
I(x_i+Δx,y_i+Δy)=I(x_i,y_i )+I_x (x_i,y_i )  Δx+I_y (x_i,y_i )  Δy
I(x_i+Δx,y_i+Δy)=I(x_i,y_i )+[I_x (x_i,y_i )     I_y (x_i,y_i )][■(Δx@Δy)]

	Remplaçons maintenant dans la fonction d’autocorrélation :
c(x,y)=∑_W▒[I(x_i,y_i )-I(x_i,y_i )-[I_x (x_i,y_i )     I_y (x_i,y_i )][■(Δx@Δy)]]^2 
	Ceci on peut le simplifier :
c(x,y)=[Δx  Δy][■(∑_W▒I_x^2 &∑_W▒〖I_x I_y 〗@∑_W▒〖I_x I_y 〗&∑_W▒I_y^2 )][■(Δx@Δy)]
	La matrice :
M=[■(∑_W▒I_x^2 &∑_W▒〖I_x I_y 〗@∑_W▒〖I_x I_y 〗&∑_W▒I_y^2 )]
	Illustre la structure du gradient de l’image dans un voisinage du pixel (x,y). 
	Ses valeurs propres et décrivent l’étalement des valeurs du gradient dans deux directions orthogonales 	(et peuvent être associées aux “courbures” principales de la fonction d’autocorrélation).

Etape 3 : On calcule la valeur R (la réponse de l’opérateur Harris) définie comme :
R=det⁡(M)-k×trace(M)^2
R=λ_1 λ_2-k(λ_1+λ_2 )^2
	k : le facteur de sensibilité pour séparer les coins des bords, typiquement une valeur proche 	de zéro 	(0.4 ≤ k ≤ 0.6)
Etape 4 : L’analyse des iso-contours de R fournit l’interprétation suivante : 
 
Méthode de Förstner
Förstner identifie également les points d’intérêt à l'aide de la matrice M.
Or, l’algorithme de Förstner est utilisé pour obtenir une solution approximative avec une précision sous-pixel. Il résout le point le plus proche de toutes les lignes tangentes du coin dans une fenêtre donnée et est une solution par moindres carrées. L'algorithme repose sur le fait que pour un coin idéal, les lignes tangentes se croisent en un seul point.
Etape 1 : Calcul des images du gradient I_x et I_y
Les filtres de Roberts sont une approche discrète de la dérivée de pas 1 d'une fonction : le gradient de cette fonction.
Si I(x,y) represente un pixel dans une image, alors les amplitudes des gradients en x et en y peuvent s'ecrire respectivement :
I_x=∂I/∂x= I(x_(i+1),y_(j+1)) - I(x,y)
I_y=∂I/∂y= I(x,y_(j+1)) - I(x_(i+1),y)
Cela revient à convoluer l'image avec les deux filtres R_x = [■(+1&0@0&-1)] et R_y = [■(0&+1@-1&0)]
Etape 2 : pour tous les points de l’image de gradient, calculer la matrice de covariance du gradient, et en extraire les valeurs propres λ_1 et λ_2
M=[■(∑_W▒I_x^2 &∑_W▒〖I_x I_y 〗@∑_W▒〖I_x I_y 〗&∑_W▒I_y^2 )]
Etape 3 : Calcul des valeurs d’intérêts q et w avec M
Contrairement à Harris, Förstner prend en compte les deux valeurs propres λ_1 et λ_2 de l'inverse de la matrice comme valeur d'intérêt. Elles définissent les axes d'une ellipse d'erreur. Par le calcul de leur taille :
w=(λ_1 λ_2)/(λ_1+λ_2 )=(det⁡(M))/(trace(M));     w>0
Et le facteur de forme relatif à la circularité :
q=1-((λ_1-λ_2)/(λ_1+λ_2 ))^2=(4 det⁡(M))/trace(M)²;      0≤q≤1
Les propriétés suivantes peuvent être déduites :
	Les petites ellipses circulaires définissent un point d’intérêt
	Les ellipses d'erreur allongées suggèrent un bord droit
	Les grandes ellipses marquent une zone homogène
Etape 4 : Détermination des points candidats avec la condition q>T_q et w>T_w
	Un point d'intérêt est présent si les valeurs seuils données T_w et T_q sont dépassées. Les paramètres 	appropriés pour cela se situent dans l'intervalle :
{█(T_w=(0.5…1.5)  w ̅@T_q=(0.5…0.75))┤
	Avec w ̅ la moyenne de toutes les valeurs w de l’image.
Etape 5 : Déterminer les coins comme points extrêmes qui prennent les plus grandes valeurs des points candidats
