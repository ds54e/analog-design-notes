"""MIT: finite local HTML/resource link check under the Pages project prefix."""
import argparse,json
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote,urlsplit

class Links(HTMLParser):
    def __init__(self):super().__init__();self.urls=[]
    def handle_starttag(self,tag,attrs):
        for key,value in attrs:
            if key in ('href','src') and value:self.urls.append(value)

def main():
    p=argparse.ArgumentParser();p.add_argument('--site',type=Path,default=Path('site/_site'));a=p.parse_args();root=a.site.resolve();broken=[];count=0
    for file in root.rglob('*.html'):
        parser=Links();parser.feed(file.read_text())
        for url in parser.urls:
            u=urlsplit(url)
            if u.scheme or u.netloc or not u.path:continue
            path=unquote(u.path)
            if path.startswith('/analog-design-notes/'):target=root/path[len('/analog-design-notes/'):]
            elif path.startswith('/'):broken.append((str(file.relative_to(root)),url,'outside project prefix'));continue
            else:target=file.parent/path
            target=target.resolve();count+=1
            if not target.is_relative_to(root) or not target.exists():broken.append((str(file.relative_to(root)),url,'missing'))
    required=['index.html','methods.html','ja/gain-operating-point-accuracy.html','studies/gain-operating-point-accuracy/index.html','studies/gain-operating-point-accuracy/reproduce.html','downloads/gain-study-v1.0.zip','revision.json']
    for name in required:
        if not (root/name).is_file():broken.append(('required',name,'missing'))
    assert not broken,broken
    print(json.dumps(dict(status='PASS',local_links_checked=count,required_files=len(required))))

if __name__=='__main__':main()
