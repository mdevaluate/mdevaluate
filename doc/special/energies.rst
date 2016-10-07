Gromacs Energy Files
====================

It is possible to read the energy files (.edr) written out by Gromacs with mdevaluate.
Those files contain thermodynamic properties of the system, like temperature or pressure.
The exact contents of an energy file depend on the type of ensemble that was simulated,
an NVT simulation's energy file for example will not contain information about the box size.

To open these files use the function :func:`mdevaluate.open_energy`, which takes the filename of an energy file.
The types of energies stored in the file can be shown with the :attr:`types` attribute of the class :class:`~mdevaluate.reader.EnergyReader`,
the :attr:`units` attribute gives the units of these energy types.
The timesteps at which those energies were written out are accessible through the :attr:`~mdevaluate.reader.EnergyReader.time` property.
The time series of one of these energies can be accessed through the named index, comparable to python dictionaries.
::
  import mdevaluate as md
  edr = md.open_energy('/path/to/energy.edr')
  # plot the evolution of temperature
  plot(edr.time, edr['Temperature'])
