"""MIT: finite local HTML/resource link check under the Pages project prefix."""
import argparse,json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote,urlsplit

class Links(HTMLParser):
    def __init__(self):super().__init__();self.urls=[];self.ids=set()
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ('href','src') and value:self.urls.append(value)
            if key=='id' and value:self.ids.add(value)

def main():
    p=argparse.ArgumentParser();p.add_argument('--site',type=Path,default=Path('site/_site'));a=p.parse_args();root=a.site.resolve();broken=[];count=0
    parsed={}
    for file in root.rglob('*.html'):
        parser=Links();parser.feed(file.read_text());parsed[file]=parser
    for file,parser in parsed.items():
        for url in parser.urls:
            u=urlsplit(url)
            if u.scheme or u.netloc:continue
            if not u.path:
                if u.fragment and unquote(u.fragment) not in parser.ids:broken.append((str(file.relative_to(root)),url,'missing local anchor'))
                continue
            path=unquote(u.path)
            if path.startswith('/analog-design-notes/'):target=root/path[len('/analog-design-notes/'):]
            elif path.startswith('/'):broken.append((str(file.relative_to(root)),url,'outside project prefix'));continue
            else:target=file.parent/path
            target=target.resolve();count+=1
            if not target.is_relative_to(root) or not target.exists():broken.append((str(file.relative_to(root)),url,'missing'))
            elif u.fragment and target.suffix=='.html' and unquote(u.fragment) not in parsed[target].ids:broken.append((str(file.relative_to(root)),url,'missing target anchor'))
    required=['index.html','methods.html','ja/gain-operating-point-accuracy.html','studies/gain-operating-point-accuracy/index.html','studies/gain-operating-point-accuracy/reproduce.html','downloads/gain-study-v1.0.zip','revision.json']
    required += ['ja/passive-sensitivity.html','studies/passive-sensitivity/index.html','studies/passive-sensitivity/reproduce.html','downloads/passive-study-v1.0.zip']
    required += ['ja/active-allocation.html','studies/active-allocation/index.html','studies/active-allocation/reproduce.html','downloads/allocation-study-v1.0.zip']
    required += ['ja/finite-feedback.html','studies/finite-feedback/index.html','studies/finite-feedback/reproduce.html','downloads/feedback-study-v1.0.zip']
    required += ['ja/power-sequencing.html','studies/power-sequencing/index.html','studies/power-sequencing/reproduce.html','downloads/sequencing-study-v1.0.zip']
    required += ['learn/index.html','learn/01-mos-operating-point.html','learn/02-common-source.html',
        'learn/03-mirrors-and-bias.html','learn/04-differential-pair-and-ota.html','learn/05-feedback-and-response.html',
        'learn/06-errors-and-resources.html','learn/07-supply-and-startup.html','learn/08-toward-an-ldo.html',
        'ja/learning-path.html','studies/learning-path-v1/reproduce.html','downloads/learning-path-v1.0.zip']
    required += ['studies/learning-day-v1/reproduce.html','downloads/learning-day-v1.0.zip',
                 'learn/figures/initial-d1-load-response.svg','licenses/LICENSES/APM-Apache-2.0.txt']
    for name in required:
        if not (root/name).is_file():broken.append(('required',name,'missing'))
    assert not broken,broken
    print(json.dumps(dict(status='PASS',local_links_checked=count,required_files=len(required))))

if __name__=='__main__':main()
