# Project structure

## Intentions

Main intentions of the project are to:

* learn about working on bioinformatic topics,
* implement solutions in some organized way for simpler reusability, or at least to easier find what was done in the first place,
* have functionalities that are possible to mix, to allow working on problem variations,
* have possibility to run (demanding) calculation locally, on dedicated server, or on cluster.


## Command

Tool is implemented through a script `bin/zcit.py`. Script general format is
```
zcit.py <command> <arguments>
```

General help is printed with:
```
zcit.py help
```

Help for specific command is printed with:
```
zcit.py help <command>
```


## Structure

Idea is that a research is done in steps. Steps have similar data format and functionalities between researchs.
Step depends on previous step(s) and/or human decisions.

### Research

Research (project) is a directory. It is created with command `zcit.py init dirname`

### Steps

Step is a subdirectory of a research. Step has data type and it's stored data.
Step's data type is decided when step is created.
Step stored data and possible processing of that data depends on step's data type.

Example of data types: table, set of sequences, set of sequences with annotations, alignment, ...

Note: step is implemented as [abstract data type](https://en.wikipedia.org/wiki/Abstract_data_type).

Step stores 'history information' how it is created, to track research status.

### Workflow

Workflow is a research with defined step structure. It is created with `init` command by providing workflow name and parameters.
It is used to simplify work on researches with same pipeline of steps.

When research is created as a workflow, than command `zcit.py workflow` can be used to run od inquire research.
```
# Run whole research
zcit.py workflow run

# Show graph of research status
zcit.py workflow graph
# Note: it is also possible to use standard graph presentation
zcit.py graph

# Create research summary
zcit.py workflow summary

```


## Code structure

ToDo: ...
