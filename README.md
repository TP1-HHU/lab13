# Lab 11
### Nonlinear Schrödinger equation
----
Implement a splitting method for the nonlinear Schrödinger equation

<p align="center">
<img src="stuffy_stuff/f1.png" width="250">
</p>

Use an implicit scheme for the step of the linear system
<p align="center">
<img src="stuffy_stuff/f2.png" width="170">
</p>
and make explicit use of the conservation of
<img src="stuffy_stuff/f3.png" width="25"> in the sub-step for the nonlinear system
<p align="center">
<img src="stuffy_stuff/f4.png" width="170">
</p>

---
Plot in GNUplot using
`do for [i=0:50]{p 'psi_'.i; pause 0.1;set yr[0:0.4]} `
