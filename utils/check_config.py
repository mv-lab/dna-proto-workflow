import yaml
import yamlordereddictloader

def readconfig (myfile):
    with open(myfile) as f:
        yaml_data = yaml.load(f, Loader=yamlordereddictloader.Loader)
        return dict(yaml_data)

def pconfig (myconfig):
    try:
        for key, value in myconfig.items():
            print ('\n> ',key,'\n',value)
    except AttributeError:
        print ('Invalid config. NOT dictionary')

