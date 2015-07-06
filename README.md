Compiler

gcc-5 -fopenmp -O3 viho_alt.c -lfftw3 -lX11 -ltiff -ljpeg -lpng

compile with openMP



dans la version alternative, non interactive, on donne l'image puis les 9 coeffs ; viho_alt sauvegarde le résultat dans une nouvelle image (au format .png, de taille WOUTxHOUT = 512x512)

dans la version démo, non interactive, on donne l'image, puis quatre points (au format abscisse ordonnée) puis les quatre points sur lesquels ont souhaite les envoyer ; viho_demo sauvegarde le résultat dans une nouvelle image (au format .png, de taille WOUTxHOUT = 512x512)
exemple : tissu_a est un tapis dont les sommets, listés dans le sens horaire, ont pour coordonnées (495,5), (898,77), (596,634) et (10,261) donc l'appel
  ./a tissu_a.jpg 495 5 898 77 569 634 10 261 0 0 512 0 512 512 0 512
  redresse le tapis en un carré 512x512

dans la version warp to rectangle, non interactive, on donne l'image, puis quatre points définissant un quadrilatère convexe (sinon message d'erreur) et un ratio; viho warp to rectangle réalise l'homographie "la plus faible" (comprendre "qui tourne le moins les diagonales") transformant ce quadrilatère en rectangle, puis sauvegarde le résultat dans une nouvelle image (au format .png). La taille du rectangle de sortie est tel que hOut = ratio x wOut, et les aires du quadrilatère de départ et du rectangle d'arrivée sont égales.
exemple : 
  ./a tissu_a.jpg 495 5 898 77 569 634 10 261 1.
  et
  ./a tissu_a.jpg 898 77 10 261 495 5 569 634 1.
  redresse le tapis en un carré de même aire que le tapis
