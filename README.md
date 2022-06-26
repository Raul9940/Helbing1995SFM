# Helbing1995SFM
Codes to run simulations under different configurations of the original Social Force Model introduced in 1995 by D.Helbing and P.Molnár [1].

The following simulations are available:

-[Empty corridor](<Codes/empty_corridor.py>): Pedestrians move from one end of the corridor to the other in counterflow.

-[Width decrease in the middle of the corridor](<Codes/middle_aperture.py>): A reduction in the width at the central point of the corridor takes place. Furthermore, each group of pedestrians (moving from left to right and moving from right to left) are made to cross the aperture through a different area. They obey with a given probabilty.

-[Two reduced apertures in the middle of the corridor](<Codes/middle_aperture_with_central_wall.py>): Two reduced apertures are located in the middle of the corridor. Each group of pedestrians is made to cross through a different one. They obey with a given probability.

Simulations are run in an Euler scheme and no periodic boundary conditions are used.

References:

[1] Dirk Helbing and Péter Molnár. Social force model for pedestrians dynamics. Phys. Rev. E, 51:4282-4286. May 1995.
