import quantum.enzyme_kinetics as ek, inspect, textwrap
print('enzyme_kinetics.py =', inspect.getsourcefile(ek))
print('has _to_float?    =', hasattr(ek, '_to_float'))
print('uses odeint?      =', 'odeint' in inspect.getsource(ek))
# Print the first lines of integrate() to verify numeric casts
src = inspect.getsource(ek.EnzymeKinetics.integrate)
print(textwrap.dedent('\n'.join(src.splitlines()[:30])))
