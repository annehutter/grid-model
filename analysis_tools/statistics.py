import numpy as np

def calc_bubble_size_distribution(field, gridsize, boxsize, threshold, Nrays):
    
    assert(np.shape(field) == (gridsize, gridsize, gridsize))
    
    pos_threshold = np.where(field > threshold)
    numCells_threshold = len(pos_threshold[0]) 
    
    if(numCells_threshold == 0):
        distance = [ None for y in range(gridsize) ]
        histogram = [ None for y in range(gridsize) ]
        return distance, histogram
    else:
        histogram = np.zeros(gridsize)
        pos = np.zeros(3, dtype=np.int32)
        radii = []
            
        for ray in range(Nrays):
            index = np.random.randint(low=0, high=numCells_threshold)
            
            # pick point
            for i in range(len(pos)):
                pos[i] = pos_threshold[i][index]
            
            if(field[pos[0]][pos[1]][pos[2]] < threshold):
                print field[pos[0]][pos[1]][pos[2]]
            assert(field[pos[0]][pos[1]][pos[2]] >= threshold)
            
            # pick random direction
            direction = np.random.random_integers(low=0, high=2)
            posneg = 2*np.random.random_integers(low=0, high=1)-1
            
            distance = 0
            threshold_tmp = field[pos[0]][pos[1]][pos[2]] 
            while(threshold_tmp >= threshold and distance < gridsize-1):
                distance = distance + 1
                pos[direction] = (pos[direction] + posneg) % gridsize
                threshold_tmp = field[pos[0]][pos[1]][pos[2]]

            histogram[distance] = histogram[distance] + 1
            radii.append(np.float64(distance))

        NB = np.arange(gridsize) + 0.5
        
        factor = boxsize / gridsize
        radii = np.float64(radii) * factor
        NB = NB * factor
        
        (histogram, binedges) = np.histogram(radii, bins=NB, density = True)
        distance = 0.5 * (binedges[:-1] + binedges[1:])

        return distance, histogram
