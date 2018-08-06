FUNKTIONSUMFANG
- Matrixoperationen: Addition, Subtraktion, skalare Multiplikation,
  Matrix/Matrix Multiplikation, Determinante, ...
- Lösen von linearen Gleichungssystemen mittels QR-Zerlegung
- Teilen eines Bildes in zwei Flächen durch ein Polynom vom Grad n
  (Dabei wird auf die Distanz der Durchschnittsfarbe zwischen geteiltem
  Bild und dem Ursprungsbild optimiert.)
- Teilen eines Bildes in Quadtrees (vier Unterrechtecke) anhand eines
  Parameters für den Wert der Distanz zur Durchschnittsfarbe.



BESCHREIBUNG DER DATEIEN

CImg.h
Grafikbibliothek zum Laden und Schreiben von .bmp Bildern.
(Quelle: http://cimg.eu)

Color.h / Color.cpp 
Klasse für Farboperationen, z.b. wichtig die Distanzfunktion:
Ist nicht einfach Vektor minus Vektor, sondern enthält eine 
Annhäherung an das Farbempfinden unseres Auges.

ColorMatrix.h / ColorMatrix.cpp
Erbt von Matrix<Color>, plus die Funktionen für das Quadtree-Verfahren
sowie die Approximation durch ein Polynom.

fastbuildnrun
Hatte ich zum Entwickeln benutzt (will ja nicht immer alles kompilieren müssen)

fullbuild
Kompiliert und Link alle Files

IOcolorMatrix.cpp / IOcolorMatrix.h
Zwei Hilfsfunktionen welche ein .bmp direkt als Matrix<Color> lesen und schreiben.

main
Die Kompilierte Applikation (für macOS 10.13.4)

main.cpp
Parameterhandling und Aufruf der eigentlichen Funktionalität in ColorMatrix.

Matrix.cpp / Matrix.h
Matrix Klasse mit Template. Enthält alle gängigen Matrixoperationen wie Addition,
Multiplikation (mit Skalar und Matrix/Matrix), sowie eine Funktion zum Berechnen
der Determinante sowie weiteres.

MatrixGenerator.cpp / MatrixGenerator.h
Statische Funktionen welche bestimmte Matrizen generieren, unteranderem:
Einheitsmatrix, Givens-Rotation, Zeiletauschen, etc...

png2bmp
Kleines Shell-Script welches mittels ImageMagick PNG's in BMP's umwandelt.

Solver.cpp / Solver.h
Statische Funktion zum Lösen von linearen Gleichungssystem mittels QR-Zerlegung.



GEBRAUCHSANWEISUNG
Wichtig ist, dass jeweils nur quadratische Bilder mit einer 2er Potenz als Kantenlänge
verwendet werden (ansonsten klappt das mit den Quadtrees schlecht).
Auch können nur .bmp Bilder ohne Transparenz verwendet werden (hat mit der externen
Grafikbibliothek und dem Installationsort von ImageMagick zu tun, wenn die Bibliotheke
auf ImageMagick zugreifen könnte, so könnten auch alle Arten von Bilder geladen werden.)
Die Reihenfolge der Paramter ist wie folgt:
1. Pfad für den Bildinput
2. Pfad für den Bildoutput
3. Grenzewert für die Distanz der Durchschnittsfarbe (ist der Wert für den geteilten
   Bildausschnitt größer als dieser Wert, so wird weiter unterteilt)
4. Grad des Polynoms (-1: falls kein Polynom verwendet werden soll, es wird dann bloß
   der Durchschnittsfarbwert gemalt)
5. Maximale Größe der Bildunterteilung (solange größer als dieser Wert so wird automatisch
   ohne zu rechnen weiter unterteilt)
6. Minimale Größe der Bildunterteilung (es werden keine kleinere Unterteilungen als dieser
   Wert zulässt gemacht)



BEISPIELE
Einige Beispiele und deren Parameter zum Ausführen:

- Eiger, Mönch und Jungfrau als Polynom vom Grad 10:
  ./main input/eigermoenchjungfrau.bmp output/eigermoenchjungfrau.bmp 1000 10 512 512

- Quadtrees only:
  ./main input/sonnenuntergang.bmp output/sonnenuntergang.bmp 0.2 -1 256 4

- Ein parabolisches Blätterdach:
  ./main input/baum.bmp output/baum.bmp

- Nur mit Geraden:
  ./main input/pamplona.bmp output/pamplona.bmp 0.4 1 128 8




REPOSITORY LINK
https://github.com/ruediluethi/Polynomial_Quadtrees