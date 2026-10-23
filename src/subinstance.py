"""Bounded exact subinstance shortcut, followed by full-input verification.

A complete list for a public subset contains every full-input solution.
A censored subset search is never returned as a complete answer. The fixed
node budget is applied by the experiment runner.
"""
from hnp_native import decode,key,prepare,recover_pruned,unique
from ec_native import curve,point


def public_subset(public):
    limit={'mihnp':12,'echnp':8,'hnphm':5}.get(public['family'],public['N'])
    n=min(public['N'],limit)
    subset=dict(public);subset['N']=n
    if public['family']=='hnphm':
        fulln=public['N']
        subset['reports']=public['reports'][:n+1]+public['reports'][fulln+1:fulln+n+1]
    else:
        subset['reports']=public['reports'][:n]
        field='shifts' if public['family']=='mihnp' else 'queries'
        subset[field]=public[field][:n]
    return subset


def verify_candidates(public, candidates):
    answers=[]
    for candidate in candidates:
        if public['family']=='mihnp':
            z=(candidate+public['shifts'][0])%public['p']
            if not z:continue
            a=pow(z,-1,public['p'])-public['reports'][0]
        elif public['family']=='echnp':
            q=point(curve(public),candidate)+point(curve(public),public['queries'][0])
            if q.is_zero():continue
            a=int(q[0])-public['reports'][0]
        else:a=candidate-public['reports'][0]
        answers.extend(decode(public,a))
    return unique(answers)
