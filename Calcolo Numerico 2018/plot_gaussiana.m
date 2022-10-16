clear
close all
clc
%matlaba diegna le funzioni per punti quindi fissiamo np numero di punti
np = 100;
x=linspace(-5,5,np);
%linspace non specifico il passo, creo un vettore di punti equispaziati tra
%i primi due termini, e questo intervallo viene diviso in np spazi della
%stessa distanza che non ci interessa
%y=exp(-x.^2); - gaussiana
y=x.*sin(1./x);
plot(x,y)
%se aumento il numero di punti per cui deve passare il grafico allora devo
%aumentare np