# Usage

## Using Iggyptop

The iggytop package is not yet published so if it was not installed automatically when running 'uv sync' you can install it with 

`pip install .`

to get started right away, try running 

`uv run create_knowledge_graph.py`

which will create a knowledge graph with all available databases, convert it to AIRR format and save it to a json file.

To use just a subset of the available databases, you can just comment out adapters in the "adapters" list (in create_knowledge_graph.py)


## Additional Examples
Please check out the tutorials for more use cases.