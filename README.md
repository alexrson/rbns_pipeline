# RNA Bind-n-seq pipeline


RNA Bind-n-Seq is a method for measuring the binding affinity of proteins for specific RNA sequences. 

The paper is here: http://www.cell.com/molecular-cell/abstract/S1097-2765%2814%2900327-X 

Installation
============
To install the pipeline, do something like the following.

Get the source code
```
git clone https://github.com/alexrson/rbns_pipeline
cd rbns_pipeline
```

Install pip, if necessary. 
See instructions at http://pip.readthedocs.org/en/latest/installing.html

Install virtualenv
```
pip install virtualenv
```

Create a virtual environment for the pipeline
```
virtualenv venv 
source venv/bin/activate
```

Install the pipeline's dependencies
```
pip install -r requirements.txt
```


