package scala.auth

import pdi.jwt.{Jwt, JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._
import org.scalatestplus.play._
import auth.Auth
import auth.RSA._
import java.security._

import models.User

import scala.setup.{WithDatabase, WithToken}



class AuthSpec extends PlaySpec with WithToken with WithDatabase with OneAppPerSuite {

  /* To generate id_rsa, id_rsa.pub test files: */
  //writetestSample1Pair(512)

  val db = dbContext.db
  val Auth = new Auth(db)
  val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")

  "JWT decoding" should {

    // A few tests, but let's try pauldijou's library on this

    "return the claim when the JWT is valid" in {
      val token = Jwt.encode("""{"user":1}""", privateKey, JwtAlgorithm.RS256)
      JwtJson.decodeJson(token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "decode an auth0 token" in {
      JwtJson.decodeJson(auth0Token, auth0PublicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "fail with the wrong public key" in {
      JwtJson.decodeJson(auth0Token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Failure[_]]
    }

  }

  "validateToken" should {

    "decode an RSA token and return the claim as a User instance" in {
      val claim = Auth.validate(auth0Token, db)
      claim mustBe a[Success[_]]
      claim.get mustBe a[User]
    }

    "decode a HMAC token and return the claim as a User instance" in {
      val claim = Auth.validate(hmacToken, db)
      claim mustBe a[Success[_]]
      claim.get mustBe a[User]    }

    "fail if the token was modified, so that the signature does not match" in {
      // "name" changed
      val wrongToken = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJlbWFpbCI6Imp1bGllbi5kZWxhZm9udGFpbmVAeWFuZGV4LmNvbSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJ1c2VyX2lkIjoiYXV0aDB8NTg5NDhlMzViMGU4YzAwZjc0YjM4OTdkIiwiY2xpZW50SUQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsInBpY3R1cmUiOiJodHRwczovL3MuZ3JhdmF0YXIuY29tL2F2YXRhci9iOGJmYTc4NjY1ZGM2NzBhN2NiZjM5ODI0MzJiMmJiOD9zPTQ4MCZyPXBnJmQ9aHR0cHMlM0ElMkYlMkZjZG4uYXV0aDAuY29tJTJGYXZhdGFycyUyRmp1LnBuZyIsIm5pY2tuYW1lIjoianVsaWVuLmRlbGFmb250YWluZSIsImlkZW50aXRpZXMiOlt7InVzZXJfaWQiOiI1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJwcm92aWRlciI6ImF1dGgwIiwiY29ubmVjdGlvbiI6IlVzZXJuYW1lLVBhc3N3b3JkLUF1dGhlbnRpY2F0aW9uIiwiaXNTb2NpYWwiOmZhbHNlfV0sInVwZGF0ZWRfYXQiOiIyMDE3LTAzLTMwVDA4OjQ5OjAxLjY2N1oiLCJjcmVhdGVkX2F0IjoiMjAxNy0wMi0wM1QxNDowNTo0MS4zMzlaIiwibGFzdF9wYXNzd29yZF9yZXNldCI6IjIwMTctMDMtMjlUMTA6Mjk6NDQuMDcyWiIsIm5hbWUiOiJwaXJhdGVAZ21haWwuY29tIiwiaXNzIjoiaHR0cHM6Ly9qZGVsYWZvbi5ldS5hdXRoMC5jb20vIiwic3ViIjoiYXV0aDB8NTg5NDhlMzViMGU4YzAwZjc0YjM4OTdkIiwiYXVkIjoiRzRiRlducGF1bFJsSFNjeXp6bnhOWmVQak54WkcyOHkiLCJleHAiOjMxNDkwODYzNzQxLCJpYXQiOjE0OTA4NjM3NDEsIm5vbmNlIjoidHR0dCJ9.dc_jJDNRA7A4B6gy6jzs10t5Dn7wpF9ebU_noC5YWDt9h4HcT_W2n8NJ5jdffDXpDgFcjtkkw7NDEj5OlLz6pAPdBL8cqAqdxXWWHXCxomcLerOpAFjxXCikBNnAntQklfxnazJ28OEXzqzZKuSjoGO_kECewubhT-j-8C90wAU"
      Auth.validate(wrongToken, db) mustBe a[Failure[_]]

      // "iss" changed
      val wrongToken2 = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJlbWFpbCI6Imp1bGllbi5kZWxhZm9udGFpbmVAeWFuZGV4LmNvbSIsImVtYWlsX3ZlcmlmaWVkIjp0cnVlLCJ1c2VyX2lkIjoiYXV0aDB8NTg5NDhlMzViMGU4YzAwZjc0YjM4OTdkIiwiY2xpZW50SUQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsInBpY3R1cmUiOiJodHRwczovL3MuZ3JhdmF0YXIuY29tL2F2YXRhci9iOGJmYTc4NjY1ZGM2NzBhN2NiZjM5ODI0MzJiMmJiOD9zPTQ4MCZyPXBnJmQ9aHR0cHMlM0ElMkYlMkZjZG4uYXV0aDAuY29tJTJGYXZhdGFycyUyRmp1LnBuZyIsIm5pY2tuYW1lIjoianVsaWVuLmRlbGFmb250YWluZSIsImlkZW50aXRpZXMiOlt7InVzZXJfaWQiOiI1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJwcm92aWRlciI6ImF1dGgwIiwiY29ubmVjdGlvbiI6IlVzZXJuYW1lLVBhc3N3b3JkLUF1dGhlbnRpY2F0aW9uIiwiaXNTb2NpYWwiOmZhbHNlfV0sInVwZGF0ZWRfYXQiOiIyMDE3LTAzLTMwVDA4OjQ5OjAxLjY2N1oiLCJjcmVhdGVkX2F0IjoiMjAxNy0wMi0wM1QxNDowNTo0MS4zMzlaIiwibGFzdF9wYXNzd29yZF9yZXNldCI6IjIwMTctMDMtMjlUMTA6Mjk6NDQuMDcyWiIsIm5hbWUiOiJqdWxpZW4uZGVsYWZvbnRhaW5lQHlhbmRleC5jb20iLCJpc3MiOiJwaXJhdGUiLCJzdWIiOiJhdXRoMHw1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJhdWQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsImV4cCI6MzE0OTA4NjM3NDEsImlhdCI6MTQ5MDg2Mzc0MSwibm9uY2UiOiJ0dHR0In0.v8zFYyHsWDttZJLHBJJo4QZX5116HKDFBJgTc8EBfw0T3NyOwHlN6UtuagOSYQmNjKk_xMD1pgptihRcqhRTN6p5VLjCOyrga1vpxWTUfKIzbVZ2XtwHzQWhyN2qb7ipLBhXxqTJgrOgeoYA9Tb2lvDIwPW1Qq4kVaeYtkarrX8"
      Auth.validate(wrongToken2, db) mustBe a[Failure[_]]
    }

  }

  "Auth routes" should {

    "return 401 and the error message if authorization failed" in {
        val header = ("Authorization", "Bearer " + auth0Token.slice(0, auth0Token.length-1))
        val response = route(app, FakeRequest(GET, "/bai/ttt").withHeaders(header)).get
        status(response) mustBe UNAUTHORIZED
        contentAsString(response) must startWith("Last unit does not have enough valid bits")
    }

  }

}
