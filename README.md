Compiler

gcc-5 -fopenmp -O3 viho_alt.c -lfftw3 -lX11 -ltiff -ljpeg -lpng

compile with openMP



dans la version alternative, non interactive, on donne l'image puis les 9 coeffs ; viho_alt sauvegarde le résultat dans une nouvelle image (au format .png, de taille WOUTxHOUT = 512x512)

dans la version démo, non interactive, on donne l'image, puis quatre points (au format abscisse ordonnée) puis les quatre points sur lesquels on souhaite les envoyer ; viho_demo sauvegarde le résultat dans une nouvelle image

dans la version warp to rectangle, non interactive, on donne l'image, quatre points définissant un quadrilatère convexe et un ratio ; viho warp to rectangle réalise l'homographie "la plus faible" (comprendre "qui tourne le moins les diagonales") transformant ce quadrilatère en rectangle, puis sauvegarde le résultat dans une nouvelle image (au format .png). La taille du rectangle de sortie est telle que hOut = ratio x wOut, et les aires du quadrilatère de départ et du rectangle d'arrivée sont égales.
