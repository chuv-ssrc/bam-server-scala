
Authorization
=============

To ensure that only registered users have access only to the data that belongs to them,
all requests must include a Json Web Token (JWT) in their "Authorization" header.
It is a signed Bearer token which payload includes at least

* The user identifer under claim `"name"` (configurable)
* The app identifier under claim `"iss"`.

The authorization process works as follows:

1. The user ("name") is matched against the `users` table in the database
   to make sure he is registered.

2. The app name ("iss") is matched against the `apps` table in the database,
   so as to retreive the corresponding verification key from the `key` column
   and the `algorithm` used to sign the token.

3. The JWT is verified using the key.


HMAC protocol
-------------

The HMAC protocol uses a shared secret key to both encrypt and verify the encrypted data.
In our application, the authorization server encrypts the JWT signature with the secret key.

We copy the secret key in the `key` column of the `apps` table,
along with the encryption algorithm in column `algorithm`.
Bam-server can then verify token signatures using the secret key.

.. attention::

   The secret key is then as safe as your database content is.
   An attacker obtaining the key can not only forge queries for the bam-server,
   but for all other services based on the same authorization server.
   We advise to use the RSA protocol whenever possible.


RSA protocol
------------

The RSA protocol is an asymmetric algorithm that uses a private, secret key to encrypt some data,
which can be then verified (not decrypted) using a shared public key.

In our application, the authorization server encrypts the JWT signature with its private key,
and provides the public key.

We copy the public key to the bam-server database, in the `key` column of the `apps` table,
along with the encryption algorithm in column `algorithm`.
Bam-server can then verify token signatures using the public key.

RSA keys are not strings, but can be stored and shared as text in these common formats:

* PEM format:

  .. code:: bash

      -----BEGIN RSA PUBLIC KEY-----
      MFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt
      xgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==
      -----END RSA PUBLIC KEY-----

* CER format ("certificate"):

  .. code:: bash

      -----BEGIN CERTIFICATE-----
      MIIC8jCCAdqgAwIBAgIJY0fhkAkJqc1YMA0GCSqGSIb3DQEBBQUAMCAxHjAcBgNV
      ...
      -----END CERTIFICATE-----

Copy one of these into the `key` column of the `apps` table.

.. note::

  Database TEXT fields can only hold one line of text.
  You must replace all line returns by **`\\n`**, for instance:

  .. code:: bash

      -----BEGIN RSA PUBLIC KEY-----\nMFwwDQYJKoZIhvcNAQEBBQADSwAwSAJBAK7ttYaE/1ldsb0OJQDQhhDWqwuFWIyt\nxgYIJH1HYA4UpA/Nm24fERIA1xi2Pomep6VTnQ/ThFP5hn2NyITwCIsCAwEAAQ==\n-----END RSA PUBLIC KEY-----


JWT example
-----------

.. code:: bash

    # Header
    {
      "typ": "JWT",
      "alg": "RS256"
    }

    # Payload
    {
      "name": "myUsername77",
      "iss": "myAppname",
      "exp": 31490863741,
      "iat": 1490863741,
      "sub": ...,
      "aud": ...,
      ...
    }

    # Signature
    RSASHA256(
      base64UrlEncode(header) + "." +
      base64UrlEncode(payload),
      <public key>,
      <private key>
    )

    # Final token
    <Base64(Header)>.<Base64(Payload)>.<Base64(Signature)>

For more details on JWTs, see `jwt.io <jwt.io>`_.
