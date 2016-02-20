EoN

"Epidemics on Networks"
======

    EoN is a Python package for the simulation of epidemics on networks.

    We assume that networks are created using the NetworkX package.

    The algorithms are based on the book:
            Mathematics of network epidemics: from exact to approximate models
            - Kiss, Miller & Simon
    The book is not yet published, and so at the moment this is intended primarily
    as a placeholder for references from the book.  As the book nears publication, 
    a more complete readme and guide will be provided.

    Please cite the book if using these algorithms


Highlights:
---------
A highlight of this work is the use of event-based simulations: "fast_SIR" and "fast_SIS".
These will outperform Gillespie simulations, and other methods.




Further comments:
---------------
This is a preliminary version:
  - The algorithms have not been tested in Python 3 (tested only in 2.7)
  - The figure references are based on the current draft of the book, which
       may still be edited, so figure numbers may change.
  - Additional changes are still needed.   I plan to add:
     * an implementation of fast_SIS
     * A Gillespie Algorithm (SIS and SIR)
     * More complete and rigorous testing
     * A quickstart guide.
 
