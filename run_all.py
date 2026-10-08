"""Fresh-run reproduction of the manuscript's finite certificates.

Run python run_all.py --core for the two standard-library graph certificates.
Run python run_all.py for every exact check, including SymPy and NetworkX.
"""
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent

def read(p):return json.loads(p.read_text(encoding='utf-8'))

def manifest():
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        digest,name=line.split('  ',1)
        p=(ROOT/name).resolve()
        assert p.is_relative_to(ROOT),'manifest path escapes archive'
        assert hashlib.sha256(p.read_bytes()).hexdigest()==digest,'manifest mismatch: '+name

def compare_algorithms(folder,full):
    old=read(folder/'certificate_results.json')
    frontier=read(folder/'frontier_results.json')
    array=read(folder/'rank_counts.json')
    for i,w in enumerate(frontier['witnesses']):
        assert w['counts']==array
        assert w['flow']==old['witnesses'][i]['flow_coefficients']
        for x,y in zip(w['squares'],old['witnesses'][i]['squares']):
            assert x['cycle']==y['cycle']
            assert x['contracted']==y['contracted_coefficients']
            assert x['crossed']==y['crossed_coefficients']
    if full:
        dc=read(folder/'independent_audit_results.json')['flow_coefficients']
        for a,b in zip(dc,frontier['witnesses']):
            assert a['base']==b['flow']
            for x,y in zip(a['modifications'],b['squares']):
                assert x['cycle']==y['cycle']
                assert x['contracted']==y['contracted'] and x['crossed']==y['crossed']

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core',action='store_true')
    args=parser.parse_args()
    if sys.flags.optimize:raise RuntimeError('Do not run with -O: assertion checks must be enabled.')
    manifest()
    env={'python':sys.version,'platform':platform.platform(),'packages':{}}
    if not args.core:
        for package,version in [('sympy','1.14.0'),('networkx','3.6.1')]:
            actual=importlib.metadata.version(package)
            if actual!=version:raise RuntimeError(f'{package}: expected {version}, found {actual}')
            env['packages'][package]=actual
    jobs=[('verify_certificates.py',['rank_counts.json','certificate_results.json']),
          ('verify_frontier.py',['frontier_results.json'])]
    if not args.core:
        jobs += [('check_trace.py',['trace_results.json']),
                 ('verify_cubic_local.py',['cubic_local_results.json']),
                 ('verify_extensions.py',['independent_audit_results.json']),
                 ('verify_general_degree_local.py',['full_local_audit_results.json']),
                 ('verify_two_spin.py',['two_spin_results.json'])]
    output=ROOT/('run_results_core' if args.core else 'run_results')
    output.mkdir(exist_ok=True)
    result={'mode':'core' if args.core else 'full','environment':env,
            'manifest_verified':True,'jobs':[],'all_passed':False}
    child_env=os.environ.copy()
    child_env['PYTHONIOENCODING']='utf-8'
    child_env['PYTHONHASHSEED']='0'
    child_env.pop('PYTHONOPTIMIZE',None)
    try:
        with tempfile.TemporaryDirectory(prefix='scratch_',dir=output) as temp:
            scratch=Path(temp).resolve()
            assert scratch.is_relative_to(output.resolve())
            for p in ROOT.glob('*.py'):shutil.copy2(p,scratch/p.name)
            shutil.copytree(ROOT/'inputs',scratch/'inputs')
            for script,files in jobs:
                start=time.perf_counter()
                run=subprocess.run([sys.executable,script],cwd=scratch,env=child_env,
                                   stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                (output/(Path(script).stem+'.log')).write_bytes(run.stdout)
                row={'script':script,'seconds':round(time.perf_counter()-start,3),'exit_code':run.returncode}
                result['jobs'].append(row)
                if run.returncode!=0:raise RuntimeError(f'{script} failed; see its log')
                for filename in files:
                    actual=read(scratch/filename)
                    expected=read(ROOT/'expected'/filename)
                    assert actual==expected,'expected output mismatch: '+filename
                    shutil.copy2(scratch/filename,output/filename)
                row['expected_outputs_match']=True
                print('PASS',script,row['seconds'],'seconds',flush=True)
            compare_algorithms(scratch,not args.core)
            result['cross_algorithm_comparisons_passed']=True
            result['all_passed']=True
    except Exception as exc:
        result['error']=str(exc)
        raise
    finally:
        (output/'summary.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('PASS: fresh computations, expected outputs, and independent algorithms agree.')

if __name__=='__main__':main()
