Hello! This is a repo for an end-of-term Optimization project.

## Project description
The problem dealt with here was inspired by the long-running Micromouse contest (see this [Wiki](https://en.wikipedia.org/wiki/Micromouse) page for more info). It is stated like so:
Given an m*n point grid (lattice) with obstacles in the form of edges connecting any two points on the grid, find the fastest path on the grid connecting a given pair of start-target points. There are two criterion for path evaluation as follows
- Path length (in Euclidean distance)
- Number of turns made (divided into acute turns and otherwise)

(Velocity and turn costs are all assumed to be constant). For information on how the constants are decided, further comments and breakdown of ideas in solving the problem, see our presentation slides (in Vietnamese) [here]().

Here is an example of such a grid, with a start marked in blue, a target in red, and obstacles in black

![](/Sample_images/Size4/sample3.png "Example grid")

and here is a proposed solution

![](/Proposed_solutions/Visualize/Size4/sample3.png "Propsed solution")

### Brief overview of methods implemented
Two methods were attempted to solve the problem:
- The first is through a non-linear program, solved using gurobi. We were able to partially deal with the samples created, specifically the ones with feasible paths.
- The second is through graphs, solved using Dijkstra's algorithm. This is implemented on python, with help from the [dijkstar](https://pypi.org/project/Dijkstar/) library.

### Usage
The resources provided in this repo are relatively self-contained, in the sense that they can be run for existing samples. If you want to solve your own grids, save it in a .json file with content of the following format

    {
        "row": int,
        "column": int,
        "edges": list(list),
        "start": list,
        "target": list,
    }
This is a dictionary that holds the content of your grid, see the folder 'Samples' for more info. If you bother with writing and reading auxiliary files, you can tweak the codes given to generate and work with data inside the program.

### Acknowledgements
We wish to thank our teachers and classmates for their interest and assistance throughout the development of this project.

(For a Vietnamese version of the description, see the enclosed [readme](/readme(Vietnamese).md))