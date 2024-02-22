from congrid import congrid
import numpy as np

def interpshift(a, bbox_a, b, bbox_b):

    # Step 1: transpose both arrays, to handle stupid python
    # convention whereby in images, arrays are indexed as y,x rather
    # than x,y
    at=np.transpose(a)
    bt=np.transpose(b)

    # Step 2: get mesh spacing and dimensions of two inputs arrays
    sh_a=at.shape
    sh_b=bt.shape
    dx_a=np.array([(bbox_a[1]-bbox_a[0])/sh_a[0],
                   (bbox_a[3]-bbox_a[2])/sh_a[1]])
    dx_b=np.array([(bbox_b[1]-bbox_b[0])/sh_b[0],
                   (bbox_b[3]-bbox_b[2])/sh_b[1]])

    # Step 3: find the intersection of the boxes
    bbox=(max(bbox_a[0], bbox_b[0]),
          min(bbox_a[1], bbox_b[1]),
          max(bbox_a[2], bbox_b[2]),
          min(bbox_a[3], bbox_b[3]))

    # Step 4: find the indices within the input arrays
    idx_a=np.array([int(np.fix((bbox[0]-bbox_a[0])/dx_a[0])),
                    int(np.fix((bbox[1]-bbox_a[0])/dx_a[0])),
                    int(np.fix((bbox[2]-bbox_a[2])/dx_a[1])),
                    int(np.fix((bbox[3]-bbox_a[2])/dx_a[1]))])
    idx_b=np.array([int(np.fix((bbox[0]-bbox_b[0])/dx_b[0])),
                    int(np.fix((bbox[1]-bbox_b[0])/dx_b[0])),
                    int(np.fix((bbox[2]-bbox_b[2])/dx_b[1])),
                    int(np.fix((bbox[3]-bbox_b[2])/dx_b[1]))])

    # Step 5: assign new values
    for i in range(idx_b[0], idx_b[1]+1):
        imap = ((i)*dx_b[0]+bbox_b[0]-bbox_a[0])/dx_a[0]
        for j in range(idx_b[2], idx_b[3]+1):
            jmap = ((j)*dx_b[1]+bbox_b[2]-bbox_a[2])/dx_a[1]
            bt[i,j] = at[int(imap), int(jmap)]
    # Step 4: do nearest neighbor interpolation
    #print idx_a
    #print idx_b
    #b[idx_b[0]:idx_b[1]+1, idx_b[2]:idx_b[3]+1] = \
    #    congrid(a[idx_a[0]:idx_a[1]+1, idx_a[2]:idx_a[3]+1],
    #            (idx_b[1]-idx_b[0]+1, idx_b[3]-idx_b[2]+1),
    #            method='nearest', centre=False, minusone=True)

    # Step 6: undo the transpose
    b=np.transpose(bt)
