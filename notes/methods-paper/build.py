"""Build and check both complete PDFs; use only installed TeX resources.

On a normally configured TeX Live installation no configuration is changed.
If this environment has the complete TeX tree but lacks its configuration
and generated format/map, create those generated files under the local tmp
directory. No system file or user-wide TeX configuration is modified.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys

from pypdf import PdfReader

HERE = Path(__file__).resolve().parent
TEMP = HERE / 'tmp/pdfs'
OUTPUT = HERE / 'output/pdf'
VERIFY = HERE / 'verification'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_hashes():
    return {str(p.relative_to(HERE)):digest(p) for p in sorted(HERE.rglob('*'))
            if p.is_file() and p.suffix in ('.tex', '.bib')
            and not set(p.relative_to(HERE).parts) & {'tmp', 'output'}}


def locate(name, env):
    result = subprocess.run(['kpsewhich', name], env=env, capture_output=True, text=True)
    return result.stdout.strip() if result.returncode == 0 else ''


def configured_environment():
    env = dict(os.environ)
    env['SOURCE_DATE_EPOCH'] = str(int(datetime(2026, 11, 27, tzinfo=timezone.utc).timestamp()))
    configuration = {'local_format_created':False, 'local_font_map_created':False}
    if not locate('article.cls', env):
        distribution = Path('/usr/share/texlive/texmf-dist')
        extra = Path('/usr/share/texmf')
        if not (distribution / 'tex/latex/base/article.cls').is_file():
            raise RuntimeError('A complete TeX Live installation is required; article.cls is missing.')
        env.update(TEXMFCNF=str(distribution / 'web2c'), TEXMFROOT=str(distribution.parent),
                   TEXMF='{'+str(distribution)+','+str(extra)+'}')
        configuration['installed_tree_configuration'] = str(distribution)
    local = HERE / 'tmp/tex-formats'
    local.mkdir(parents=True, exist_ok=True)
    if not locate('pdflatex.fmt', env):
        env['TEXFORMATS'] = str(local)+'//:'
        if not (local / 'pdflatex.fmt').is_file():
            with (local / 'initialize.log').open('w') as log:
                subprocess.run(['pdftex','-ini','-etex','-interaction=nonstopmode',
                                '-halt-on-error','-jobname=pdflatex','-progname=pdflatex',
                                r'\pdfoutput=1 \input latex.ltx'],
                               cwd=local, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
        configuration['local_format_created'] = True
    if not locate('pdftex.map', env):
        maps = []
        for name in ('lm.map','cm.map','cmextra.map','euler.map','symbols.map'):
            path = locate(name, env)
            if not path: raise RuntimeError('Installed font map is missing: '+name)
            maps.append(Path(path).read_text())
        (local / 'pdftex.map').write_text('\n'.join(maps))
        env['TEXFONTMAPS'] = str(local)+'//:'
        configuration['local_font_map_created'] = True
    for name in ('article.cls','amsmath.sty','lmodern.sty','pdflatex.fmt','pdftex.map'):
        if not locate(name, env): raise RuntimeError('TeX configuration does not resolve '+name)
    return env, configuration


def verify_existing():
    record = json.loads((VERIFY / 'build.json').read_text())
    if record['sources_sha256'] != source_hashes():
        raise ArithmeticError('The manuscript sources changed after the recorded PDF build.')
    for name, value in record['documents'].items():
        path = HERE / value['pdf']
        if digest(path) != value['sha256'] or len(PdfReader(path).pages) != value['pages']:
            raise ArithmeticError('The delivered PDF differs from its recorded build: '+name)
        if value['overfull_boxes'] or value['unresolved_references']:
            raise ArithmeticError('The recorded build retains a layout/reference defect.')
    print(json.dumps({'verified_documents':{name:v['pages'] for name,v in record['documents'].items()},
                      'all_source_and_pdf_hashes_match':True},indent=2))


def build():
    subprocess.run([sys.executable, str(HERE / 'generate_tables.py')], check=True, cwd=HERE)
    TEMP.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    VERIFY.mkdir(parents=True, exist_ok=True)
    env, configuration = configured_environment()
    env['BIBINPUTS'] = str(HERE) + ':' + env.get('BIBINPUTS', '')
    version = subprocess.check_output(['pdflatex','--version'],env=env,text=True).splitlines()[0]
    documents = {}
    for name in ('paper','supplement'):
        command = ['pdflatex','-interaction=nonstopmode','-halt-on-error','-file-line-error',
                   '-output-directory',str(TEMP),name+'.tex']
        transcript = TEMP / (name+'_passes.log')
        with transcript.open('w') as output:
            subprocess.run(command,cwd=HERE,env=env,stdout=output,stderr=subprocess.STDOUT,check=True)
            if name == 'paper':
                subprocess.run(['bibtex',name],cwd=TEMP,env=env,
                               stdout=output,stderr=subprocess.STDOUT,check=True)
            for _ in range(3):
                subprocess.run(command,cwd=HERE,env=env,stdout=output,stderr=subprocess.STDOUT,check=True)
        log = (TEMP / (name+'.log')).read_text(errors='replace')
        overfull = re.findall(r'Overfull \\[hv]box[^\n]*',log)
        unresolved = re.findall(r'[^\n]*(?:undefined references|Citation[^\n]*undefined|Reference[^\n]*undefined|multiply defined)[^\n]*',log)
        (VERIFY / (name+'_build.log')).write_text(log)
        if overfull or unresolved:
            raise ArithmeticError(json.dumps({'document':name,'overfull':overfull,'unresolved':unresolved}))
        output = OUTPUT / (name+'.pdf')
        shutil.copyfile(TEMP / (name+'.pdf'),output)
        reader = PdfReader(output)
        if not reader.pages or any(not (page.extract_text() or '').strip() for page in reader.pages):
            raise ArithmeticError('An output page is empty or not text-extractable: '+name)
        documents[name] = {'pdf':str(output.relative_to(HERE)), 'pages':len(reader.pages),
                           'sha256':digest(output),'overfull_boxes':[], 'unresolved_references':[]}
        print(name+': '+str(len(reader.pages))+' pages, no unresolved references or overflowing boxes.',flush=True)
    record = {'schema':'methods-paper-build.v1','compiler':version,'configuration':configuration,
              'sources_sha256':source_hashes(),'documents':documents}
    (VERIFY / 'build.json').write_text(json.dumps(record,indent=2)+'\n')
    verify_existing()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',help='check the existing exact source/PDF bindings')
    args = parser.parse_args()
    verify_existing() if args.verify else build()
