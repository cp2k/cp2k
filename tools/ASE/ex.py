def test():
    from ASE.Calculators.PairPotential import PairPotential
    from ASE import ListOfAtoms, Atom
    
    a = 4.0
    d = 0.71
    h2 = ListOfAtoms([Atom('H', (0, 0, 0)),
                      Atom('H', (0, 0, d))],
                     cell=(a, a, a), periodic=True)

    h2.SetCalculator(PairPotential())
    e2 = h2.GetPotentialEnergy()
    print e2



def test2():
    import CP2KCalculator
    from ASE import ListOfAtoms, Atom
    a = 4.0
    d = 0.71
    h2 = ListOfAtoms([Atom('H', (0, 0, 0)),
                      Atom('H', (d, 0, d))],
                      cell=(a, a, a), periodic=True)

    kinds ={'KIND H':{'MASS':'2', 'ELEMENT':'H',
                      'BASIS_SET':'DZVP-GTH-PADE', 'POTENTIAL':'GTH-PADE-q1'}}
    h2.SetCalculator(CP2KCalculator.CP2KCalculator(kinds))
    

    h2.GetCartesianForces()
    
    #h2.GetCalculator().UpdateParams(kinds)
    #h2.GetCalculator().UpdateParams({})
