# Geometric-Modelling development machine

## Hogyan használd

Telepítsd a Vagrant-ot a gépedre.

Miután telepítetted, add ki a következő parancsot:

~~~~
vagrant up
~~~~

Várd ki, amíg minden lefut a parancs, majd a megnyílt VirtualBox ablakban jelentkezz be.

Felhasználónév: vagrant
Jelszó: vagrant

Bejelentkezés után add ki a következő parancsot:

~~~~
startx
~~~~

Jó munkát! :)

## Megjegyzések

Az optimlib-et használó c++ kódok fordításához ihletet meríthetsz innen:

~~~~
g++ -Wall -std=c++11 -O3 -march=native -ffp-contract=fast -I/usr/include/armadillo_bits -I/usr/local/include/optim optim_de_ex.cpp -o optim_de_ex.out -lblas -llapack -loptim
~~~~

A freeglut-ot használó c++ kódok fordításához ihletet meríthetsz innen:

~~~~
g++ main.cpp -I../freeglut/include/ -lGL -lglut -lGLU
~~~~