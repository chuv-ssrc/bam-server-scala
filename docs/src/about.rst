
What is bam-server
==================

Bam-server is a microservice with a REST API that allows to query regions of
`BAM files <https://samtools.github.io/hts-specs/SAMv1.pdf>`_
remotely and securely, according to user permissions. In particular, it

* Never moves or copies the original file, but sends only the requested bits through HTTP;
* Is compatible with any authorization server that returns `JSON Web Tokens <https://jwt.io/introduction/>`_
  (*does not work with cookies*).

A typical application is the display of local read alignments for variant calling quality control,
where the subject is human and thus the data must be kept private.

An example of client web application that uses `IGV.js <https://github.com/igvteam/igv.js/tree/master>`_
to display read alignments from protected BAM files can be found here:
`bam-server-client-react <https://github.com/jdelafon/bam-server-client-react>`_.

.. figure:: /images/demo.png
   :width: 100%
   :alt: Demo interface using bam-server to display read alignments

   Demo interface using bam-server to display read alignments.

Bam-server is written in Scala using the Play framework.

