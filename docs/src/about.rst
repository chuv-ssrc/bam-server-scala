
A simple BAM server in Scala
============================

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


Methods
-------

Bam-server in itself represents only the resource server in the authorization workflow (see figure below):

1. The client app connects to its authorization server to request an authorization token;
2. The auth server returns a JSON Web Token (JWT) signed with either RSA or HMAC.
   The JWT must contain at least an identifier for the client app (typically the *iss* claim),
   and an identifier for the user that sends the request (*name* by default, but configurable).
3. The client app requests a portion of a BAM file for a given sample by the bam-server
   using the REST API.
4. Bam-server verifies the token using either the RSA public key or the HMAC shared secret,
   then checks the user and app identifiers against its own users database.
5. Bam-server extracts the requested slice from the BAM file corresponding to the sample and
6. returns it to the client app.

.. figure:: /images/bam-server.png
   :width: 100%
   :alt: Demo interface using bam-server to display read alignments.

It is up to the auth server to define valid users. The local database of bam-server
only maps user identifiers as returned by the auth server, to sample identifiers and the corresponding BAM files.

