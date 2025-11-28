# popets-dp-for-pprl
This repository stores the data extraction and linkage scripts used in the experimental setups of the paper titled <em> The use of differential privacy for privacy-preserving record linkage: Protecting the bits but not the people </em> submitted for POPETS 2026.

This includes two sets of code:
    <ol>
        <li>The <code>data-sampler</code> directory contains the code to sample the subsets of data based on vulnerable characteristics, and generates their q-gram sets.</li>
        <li> The <code>linkage</code> directory contains the code to perform linkage on the sampled subsets, inclusive of an evaluation of the classified links.</li>
            The entry point of the linkage program is the <code>encoder.py</code> file.
    </ol>

