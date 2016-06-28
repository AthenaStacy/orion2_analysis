import h5py
import numpy as np
from optparse import OptionParser
import sys

def rebin(a, newShape):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    a = np.array(a)
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(newShape)
    evList = ['a.reshape('] + \
             ['newShape[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))

def rebinStagger(a, newShape, staggerDir):
    '''
    like rebin, but for a staggered face-centered grid
    '''
    a = np.array(a)
    ns = np.zeros(2,dtype=np.int)
    count = 0
    for i in range(0,3):
        if i==staggerDir: continue
        ns[count] = newShape[i]
        count += 1
    step = (a.shape[staggerDir]-1)/(newShape[staggerDir]-1)
    ans = np.zeros(newShape, dtype=a.dtype)
    if staggerDir == 0:
        for i in range(0, ans.shape[0]):
            a2d = a[i*step,:,:]
            a2d = rebin(a2d, ns)
            ans[i,:,:] = a2d
    elif staggerDir == 1:
        for i in range(0, ans.shape[1]):
            a2d = a[:,i*step,:]
            a2d = rebin(a2d, ns)
            ans[:,i,:] = a2d
    elif staggerDir == 2:
        for i in range(0, ans.shape[2]):
            a2d = a[:,:,i*step]
            a2d = rebin(a2d, ns)
            ans[:,:,i] = a2d
    return ans

# recursively copy the group structure
def copy_groups(key):
    for attr in fin[key].attrs.keys():
        try:
            fout.create_group(key)
        except:
            pass # the group is already created
    #
    next_keys = fin[key].keys()
    for next_key in next_keys:
        next = next_key
        if key != '/': next = key+'/'+next_key
        # if the key is a group, recursivley burrow into it.
        if type(fin[next]) == h5py.highlevel.Group: copy_groups(next)

# recursively copy all hdf5 attributes in the file hierarchy
def copy_attrs(key):
    for attr in fin[key].attrs.keys():
        dtype = h5py.h5a.open(fin[key].id,attr).dtype
        if dtype.isbuiltin: 
            fout[key].attrs[attr] = dtype.type(fin[key].attrs[attr])
        else:
            fout[key].attrs[attr] = np.array([fin[key].attrs[attr]], dtype=dtype)
    #
    next_keys = fin[key].keys()
    for next_key in next_keys:
        next = next_key
        if key != '/': next = key+'/'+next_key
        # if the key is a group, recursivley burrow into it.
        if type(fin[next]) == h5py.highlevel.Group: copy_attrs(next)

# coarsen the boxes, *:datatype=0, and *:offsets=0 datasets 
# on each level and write to fout.
def copy_datasets(coarsenRatio):
    box_idt = np.dtype([('lo_i', np.int32), ('lo_j', np.int32), ('lo_k', np.int32), 
                        ('hi_i', np.int32), ('hi_j', np.int32), ('hi_k', np.int32)])
    level_index = -1
    while(True):
        level_index += 1
        level = 'level_%d'%level_index
        try:
            keys = fin[level].keys()
        except:
            # we've run out of levels
            return
        # coarsen the cell size
        fout[level].attrs['dx'] = fin[level].attrs['dx']*coarsenRatio
        # coarsen the prob domain
        dom = fin[level].attrs['prob_domain']
        ndim = len(dom)/2
        dom = list(dom)
        for i in range(0,ndim): dom[i] = dom[i]/coarsenRatio
        for i in range(ndim,2*ndim): dom[i] = (dom[i]+1)/coarsenRatio - 1
        fout[level].attrs['prob_domain'] = np.array([tuple(dom)], dtype=box_idt)[0]
        # copy over the processor assignments
        fout[level]['Processors'] = fin[level]['Processors'].value
        # coarsen the boxes
        boxes = fin[level]['boxes'].value
        ncells = 0
        ncellsFlx = 0
        for ibox in range(0,len(boxes)):
            box = np.array(list(boxes[ibox]),dtype=np.int)
            box[0:ndim] = box[0:ndim]/coarsenRatio
            box[ndim:2*ndim] = (box[ndim:2*ndim]+1)/coarsenRatio - 1
            boxes[ibox] = tuple(box)
            # count the number of coarse cells and flux box enteries in this box
            ncells += np.product(box[ndim:2*ndim] - box[0:ndim] + 1)
            for i in range(0,ndim): 
                stagger = np.array([0,0,0],dtype=np.int)
                stagger[i] += 1
                ncellsFlx += np.product(box[ndim:2*ndim] - box[0:ndim] + 1 + stagger)
        fout[level]['boxes'] = boxes
        # loop over the various datas on this level
        for key in fin[level]:
            if not ':datatype=0' in key: continue
            if type(fin[level][key]) != h5py.highlevel.Dataset: continue
            dset_name = key.split(':')[0]
            if dset_name == 'face_data': # is there a more general way to distinguish flux boxes from cell-center data?
                # this is a flux box
                # create a dataset and offset big enough for the data
                fout[level].create_dataset(dset_name+':datatype=0',[ncellsFlx],dtype=np.float)
                fout[level].create_dataset(dset_name+':offsets=0',[len(boxes)+1],dtype=np.int64)
                # write the begining data offset
                offsets = fin[level][dset_name+':offsets=0'].value
                coarse_offset = 0
                for ibox in range(0,len(offsets)-1):
                    # read and reshape the current box data
                    box = np.array(list(boxes[ibox]),dtype=np.int)
                    # read the data
                    count = 0
                    fout[level][dset_name+':offsets=0'][ibox] = coarse_offset
                    for i in range(0,ndim):
                        # read and reshape the current data box
                        stagger = np.array([0,0,0],dtype=np.int)
                        stagger[i] += 1
                        box = np.array(list(boxes[ibox]),dtype=np.int)
                        shape = box[ndim:2*ndim] - box[0:ndim] + 1 + stagger
                        shapeFine = coarsenRatio*(box[ndim:2*ndim] - box[0:ndim] + 1) + stagger
                        n = np.product(shapeFine)
                        data = fin[level][dset_name+':datatype=0'][offsets[ibox]+count:offsets[ibox]+count+n].reshape(shapeFine,order='F')
                        count += n
                        # coarsen the data
                        data = rebinStagger(data, shape, i)
                        # flatten and write the data
                        fout[level][dset_name+':datatype=0'][coarse_offset:coarse_offset+data.size] = data.T.flatten()[:]
                        coarse_offset += data.size
                    # write the end data offset
                    fout[level][dset_name+':offsets=0'][ibox+1] = coarse_offset
            else:
                # this is a regular dataset
                # create a dataset and offset big enough for the data
                ncomps = fin[level+'/'+dset_name+'_attributes'].attrs['comps']
                fout[level].create_dataset(dset_name+':datatype=0',[ncells*ncomps],dtype=np.float)
                fout[level].create_dataset(dset_name+':offsets=0',[len(boxes)+1],dtype=np.int64)
                # read the offsets into the 1D array
                offsets = fin[level][dset_name+':offsets=0'].value
                coarse_offset = 0
                for ibox in range(0,len(offsets)-1):
                    # read and reshape the current box data
                    box = np.array(list(boxes[ibox]),dtype=np.int)
                    shape = box[ndim:2*ndim] - box[0:ndim] + 1
                    shapeFine = coarsenRatio*(box[ndim:2*ndim] - box[0:ndim] + 1)
                    shape = list(shape) + [ncomps]
                    shapeFine = list(shapeFine) + [ncomps]
                    data = fin[level][dset_name+':datatype=0'][offsets[ibox]:offsets[ibox+1]].reshape(shapeFine,order='F')
                    # coarsen the data
                    data = rebin(data, shape)
                    # flatten and write the data
                    fout[level][dset_name+':datatype=0'][coarse_offset:coarse_offset+data.size] = data.T.flatten()[:]
                    # write the data offsets
                    fout[level][dset_name+':offsets=0'][ibox] = coarse_offset
                    fout[level][dset_name+':offsets=0'][ibox+1] = coarse_offset+data.size
                    coarse_offset += data.size

###################
# input parameters, read from command line
###################
parser = OptionParser()
parser.add_option('--infile', dest='infile', 
                  help='input file name')
parser.add_option('--outfile', dest='outfile', 
                  help='output file name')
parser.add_option('--coarsenratio', dest='coarsenratio', 
                  help='The grid coarsening ratio.  Default=2', default=2)
(options, args) = parser.parse_args()
error = False
if options.infile == None: 
    error = True
    print "Error, no --infile option specified"
if options.outfile == None: 
    error = True
    print "Error, no --outfile option specified"
if error: sys.exit(-1)
###################
# begin computation
###################

fin=h5py.File(options.infile,'r')
fout=h5py.File(options.outfile,'w')
copy_groups('/')
copy_attrs('/')
copy_datasets(int(options.coarsenratio))
fin.close()
fout.close()
