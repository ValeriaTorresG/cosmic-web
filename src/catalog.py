import argparse
import numpy as np
import pandas as pd
import sys
import time as t
sys.setrecursionlimit(20000)

def web_classification(n_data, n_random, void_limit=-0.9, knot_limit=0.9):
    assert len(n_data)==len(n_random)
    r = (n_data-n_random)/(n_data+n_random)
    web = np.zeros(len(r),dtype=int)
    web[r<=void_limit] = 0
    web[(r>void_limit)&(r<=0)]    = 1
    web[(r>0)&(r<=knot_limit)]    = 2
    web[r>knot_limit]             = 3
    return web

def find_friends(first_id, all_ids, pair_ids, included_ids):
    group=[]
    loc = np.where(all_ids==first_id)[0][0]
    if included_ids[loc]: return group
    included_ids[loc]=1
    group.append(first_id)
    friends = list(pair_ids[pair_ids[:,0]==first_id,1]) + list(pair_ids[pair_ids[:,1]==first_id,0])
    for f in friends:
        group+= [f]
        group+= find_friends(f, all_ids, pair_ids, included_ids)
    group = sorted(set(group))
    return group

def find_fof_groups(pairs):
    pairs = pairs.astype(int)
    all_ids = np.unique(pairs)
    groups={}
    inc = np.zeros(len(all_ids),dtype=int)
    gid=0; total=0
    for fid in all_ids:
        fof = find_friends(fid, all_ids, pairs, inc)
        if fof:
            groups[gid]=fof
            total+=len(fof)
            gid+=1
    assert total==len(all_ids)
    return groups

def inertia_tensor(x,y,z):
    x,y,z = x-np.mean(x), y-np.mean(y), z-np.mean(z)
    r2 = x*x+y*y+z*z
    I = np.array([
      [np.sum(r2 - x*x), -np.sum(x*y),   -np.sum(x*z)],
      [-np.sum(y*x),     np.sum(r2 - y*y), -np.sum(y*z)],
      [-np.sum(z*x),     -np.sum(z*y),    np.sum(r2 - z*z)]
    ])
    vals, vecs = np.linalg.eig(I)
    idx = np.argsort(-vals)
    return np.sqrt(vals[idx]), vecs[:,idx]

def compute_group_properties(groups, pos):
    props = {k:[] for k in [
      'N','MEAN_X','MEAN_Y','MEAN_Z',
      'SIGMA_X','SIGMA_Y','SIGMA_Z','SIGMA_R',
      'LAMBDA_1','LAMBDA_2','LAMBDA_3',
      'EIGEN_1','EIGEN_2','EIGEN_3']}
    for gid, members in groups.items():
        x,y,z = pos[members].T
        # print(x)
        if len(x)>0: 
            r = np.sqrt(x*x+y*y+z*z)
            
            props['N'].append(len(members))
            props['MEAN_X'].append(x.mean()); props['MEAN_Y'].append(y.mean()); props['MEAN_Z'].append(z.mean())
            props['SIGMA_X'].append(x.std());   props['SIGMA_Y'].append(y.std());   props['SIGMA_Z'].append(z.std())
            props['SIGMA_R'].append(r.std())
            L,v = inertia_tensor(x,y,z)
            props['LAMBDA_1'].append(L[0]); props['LAMBDA_2'].append(L[1]); props['LAMBDA_3'].append(L[2])
            props['EIGEN_1'].append(v[:,0]); props['EIGEN_2'].append(v[:,1]); props['EIGEN_3'].append(v[:,2])
    return props

def make_catalog(posfile, pairfile, countfile, catalogfile, webtype,
                 void_limit=-0.9, knot_limit=0.9):
    init_t = t.time()
    pos    = np.loadtxt(posfile)
    pairs  = np.loadtxt(pairfile).astype(int)
    counts = np.loadtxt(countfile).astype(int)
    web    = web_classification(counts[:,0], counts[:,1], void_limit, knot_limit)
    if webtype=='all':
        np.savetxt(catalogfile, web, fmt='%d')
        return

    wt2id = {'void':0,'sheet':1,'filament':2,'peak':3}
    wid = wt2id[webtype]
    Np  = len(pos)
    assert len(web)==Np and pairs.min()>=0 and pairs.max()<Np

    ids = np.arange(Np)[web==wid]
    mask = np.isin(pairs[:,0],ids)&np.isin(pairs[:,1],ids)
    subp = pairs[mask]
    groups = find_fof_groups(subp)
    props  = compute_group_properties(groups, pos)
    df = pd.DataFrame.from_dict(props)
    df.to_csv(catalogfile, index=False)
    print(f'Catalog t: {init_t-t.time()} s')

if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--posfile',    required=True)
    p.add_argument('--pairfile',   required=True)
    p.add_argument('--countfile',  required=True)
    p.add_argument('--webtype',    default='all')
    p.add_argument('--catalogfile',required=True)
    p.add_argument('--void_limit',  type=float, default=-0.9)
    p.add_argument('--knot_limit',  type=float, default=0.9)
    args = p.parse_args()
    make_catalog(**vars(args))