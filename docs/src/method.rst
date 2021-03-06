
How it works
============

Bam-server in itself represents only the resource server in the authorization workflow (see figure below):

1. The client app connects to its authorization server to request an authorization token;
2. The auth server returns a JSON Web Token (JWT) signed with either RSA or HMAC.
   The JWT must contain at least an identifier for the client app (typically the *iss* claim),
   and an identifier for the user that sends the request (*name* by default, but configurable).
3. The client app requests a portion of a BAM file for a given sample from the bam-server
   using the REST API.
4. Bam-server verifies the token using either the RSA public key or the HMAC shared secret,
   then checks the user and app identifiers against its own users database.
5. Bam-server extracts the requested slice from the BAM file corresponding to the sample and
6. returns it to the client app.

.. figure:: /images/bam-server.png
   :width: 100%
   :alt: Typical authorization workflow

Notes:

* It is up to the auth server to define valid users. The local database of bam-server
  only maps user identifiers as returned by the auth server, to sample identifiers and the corresponding BAM files.

* Client apps can access bam-server remotely, but for the moment,
  bam-server itself can only access files that are located on the same machine as itself
  (the location of BAM files is defined using a configuration variable `BAM_PATH`, see :doc:`./installation`).
  One cannot provide a URL as a file path in the database.




