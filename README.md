# bact-project
Data science lab in biosciences project

## How to execute
First of all follow the instruction at [mordecai repo](https://github.com/openeventdata/mordecai) 
Install dependecies (we know: `requirements.txt` is, at the moment, garbage)
```shell
pip install -r requirements.txt
```
Run the container:
```shell
docker run -p 127.0.0.1:9200:9200 -v $(pwd)/geonames_index/:/usr/share/elasticsearch/data elasticsearch:5.5.2
```
Then execute the pipeline (it requires a lot of time):
```shell
python scrape_pathways.py
python pipeline.py
```
## TODO
- [ ] code refacoring
- [ ] documentation and code comments
- [ ] fix `requirements.txt`