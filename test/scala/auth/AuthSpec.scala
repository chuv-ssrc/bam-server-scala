package scala.auth

import pdi.jwt.{Jwt, JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import org.scalatestplus.play._
import auth.Auth._
import auth.RSA._
import java.security._
import models.User



class AuthSpec extends PlaySpec {

  /* To generate id_rsa, id_rsa.pub test files: */
  //writeTestKeyPair(512)

  val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")

  val auth0Token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJlbWFpbCI6Imp1bGllbi5kZWxhZm9udGFpbmVAeWFuZGV4LmNvbSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJ1c2VyX2lkIjoiYXV0aDB8NTg5NDhlMzViMGU4YzAwZjc0YjM4OTdkIiwiY2xpZW50SUQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsInBpY3R1cmUiOiJodHRwczovL3MuZ3JhdmF0YXIuY29tL2F2YXRhci9iOGJmYTc4NjY1ZGM2NzBhN2NiZjM5ODI0MzJiMmJiOD9zPTQ4MCZyPXBnJmQ9aHR0cHMlM0ElMkYlMkZjZG4uYXV0aDAuY29tJTJGYXZhdGFycyUyRmp1LnBuZyIsIm5pY2tuYW1lIjoianVsaWVuLmRlbGFmb250YWluZSIsImlkZW50aXRpZXMiOlt7InVzZXJfaWQiOiI1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJwcm92aWRlciI6ImF1dGgwIiwiY29ubmVjdGlvbiI6IlVzZXJuYW1lLVBhc3N3b3JkLUF1dGhlbnRpY2F0aW9uIiwiaXNTb2NpYWwiOmZhbHNlfV0sInVwZGF0ZWRfYXQiOiIyMDE3LTAzLTMwVDA4OjQ5OjAxLjY2N1oiLCJjcmVhdGVkX2F0IjoiMjAxNy0wMi0wM1QxNDowNTo0MS4zMzlaIiwibGFzdF9wYXNzd29yZF9yZXNldCI6IjIwMTctMDMtMjlUMTA6Mjk6NDQuMDcyWiIsIm5hbWUiOiJqdWxpZW4uZGVsYWZvbnRhaW5lQHlhbmRleC5jb20iLCJpc3MiOiJodHRwczovL2pkZWxhZm9uLmV1LmF1dGgwLmNvbS8iLCJzdWIiOiJhdXRoMHw1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJhdWQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsImV4cCI6MzE0OTA4NjM3NDEsImlhdCI6MTQ5MDg2Mzc0MSwibm9uY2UiOiJ0dHR0In0.fZECvmNaqCNwQXHXnzCsmsUlgvNeIIAdIkAJwFlRIaHH0E0qq9FHWVn6mSRLQLFF31gEkKA3fuQk3d7FapTmRyw8q1Fjtww4bTM9GOeTQpGnk1Q-12U5EFQ4EVcq6zPXQlspAyCjQlEQyjFvn7rxI59gRDgpG97n8hAnfLBXV0dw0yGsHcIZ_0odtHrh72ZvAIoIig6PbO5A3y9_sb4xuUQBbedltwfQ-bLSDTdJgd8gKUxfdwTE3YDLjDKn5-ekgt6heGMljGTTE7rr0qVUeMsY3fVAbKBCjsQIPkQgASRrwzatU9rRPkwZYwYy_01T8FiYR4JNs1QAiy749k22wQ"
  val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")


  "JWT decoding" should {

    "return the claim when the JWT is valid" in {
      val token = Jwt.encode("""{"user":1}""", privateKey, JwtAlgorithm.RS256)
      JwtJson.decodeJson(token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "auth0" in {
      JwtJson.decodeJson(auth0Token, auth0PublicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

  }

  "validateToken" should {

    "return the claim as a User instance" in {
      val claim = validateToken(auth0Token, "auth0")
      claim mustBe a[User]
    }

    "fail when the JWT is not valid" in {
      cancel("todo")
      //decodeJWT_RSA("what?", "secret?") mustBe a[Failure[_]]
    }

  }

}
