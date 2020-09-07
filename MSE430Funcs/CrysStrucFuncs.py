def make_supercell(cell, diff_species):
    """Append all sites in a unit cell to a structure - must have cubic lattice"""
    
    #diff_species: Boolean, if true, make different species diff colors. If false, make one basis group one color.
    #Get a copy of the structure defining the basis
    basis=cell.copy()
    superCell = cell.copy()
    
    #Create a list of all the lattice points in the cubic unit cell (i.e all the corners)
    f=[[i,j,k] for i in range(2) for j in range(2) for k in range(2)]
    
    #Remove the lattice point associated with the basis [0,0,0]
    f=f[1:]

    #Add a basis at each of the unit cell lattice points
    if diff_species:
        [superCell.append(atom.species, atom.frac_coords+f[site]) for atom in basis for site in range(len(f))]
    else:
        [superCell.append(atom.specie, atom.frac_coords+f[site]) for site in range(len(f)) for atom in basis]
    
    return(superCell)

def cubicCell(cell, a3):
    """Append all sites in a unit cell to a structure"""
    from pymatgen import Structure, Lattice
    import numpy as np
    import nglview as ngl
    
    basis=cell.copy()
    superCell = cell.copy()
    prim_vec = np.asarray(cell.lattice.as_dict().get('matrix'))
    
    #Append atoms with different combinations of lattice vectors until all coordinates exceed cubic cell [1,1,1]
    i = -1
    j = -1
    k = -1
    #Since not perfect 
    thr_a = a3*1.15
    coord_base = [0,0,0]
    new_coord = coord_base.copy()
    #while all(x <=thr_a for x in new_coord):
        #while all(x <=thr_a for x in new_coord):
            #while all(x <= thr_a for x in new_coord):
    for i in range(-3, 3):
        for j in range(-3, 3):
            for k in range(-3,3):
                new_coord = prim_vec[0]*i +prim_vec[1]*j + prim_vec[2]*k
                [superCell.append(atom.species, atom.coords + new_coord, coords_are_cartesian=True) for atom in basis]
                #k +=1
                #print(new_coord)
                #print(i, j, k)
            #j +=1
            #k=-1
            #new_coord = prim_vec[0]*i +prim_vec[1]*j + prim_vec[2]*k
        #i+=1
        #j=-1
        #k=-1
        #new_coord = prim_vec[0]*i +prim_vec[1]*j + prim_vec[2]*k

    return(superCell)

def visLattice(lattice):
    from pymatgen import Structure, Lattice
    import nglview as ngl
    unit_cell= Structure(lattice, ['Cl'], [[0,0,0]], to_unit_cell=True, coords_are_cartesian=False)
    view6=ngl.show_pymatgen(unit_cell)
    view6.clear_representations()
    view6.add_unitcell()
    return(view6)

def visUC(SC, a3):
    from pymatgen import Structure, Lattice
    import nglview as ngl
    selec=[]
    for ind, site in enumerate(SC.sites):
        if all(site.coords>=-a3*.15) & all(site.coords<=a3*1.15):
            selec=selec+[ind]
    view6 = ngl.show_pymatgen(SC)
    view6.clear_representations()
    #view6.add_representation('ball+stick', aspectRatio=10, selection=selec)
    [view6.add_representation('ball+stick', aspectRatio=5, selection=[i]) for i in selec]
    return(view6)