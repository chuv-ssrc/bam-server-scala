package scala.auth

import pdi.jwt.{Jwt, JwtAlgorithm, JwtJson}
import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import org.scalatestplus.play._
import auth.Auth._
import auth.RSA._
import java.util.Base64
import java.security._


class AuthSpec extends PlaySpec {

  /* To generate id_rsa, id_rsa.pub test files: */
  //writeTestKeyPair(512)

  val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")


  "`validateJWT`" should {

    "return the claim when the JWT is valid" in {
      val token = Jwt.encode("""{"user":1}""", privateKey, JwtAlgorithm.RS256)
      JwtJson.decodeJson(token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "auth0" in {
      val auth0Token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJpc3MiOiJodHRwczovL2pkZWxhZm9uLmV1LmF1dGgwLmNvbS8iLCJzdWIiOiJhdXRoMHw1ODk0OGUzNWIwZThjMDBmNzRiMzg5N2QiLCJhdWQiOiJHNGJGV25wYXVsUmxIU2N5enpueE5aZVBqTnhaRzI4eSIsImV4cCI6MzE0OTA3ODM4NTQsImlhdCI6MTQ5MDc4Mzg1NCwibm9uY2UiOiJ0dHR0In0.rtkxcClQjxawesOEOmau5oWE_9M1gHGSJGJ8L6KHJakFCywp6R54EatwBcL91eMtKIoN2P2G9QmVl1HB63CAI_1DwbQvBf2PzArzN6bta7GizJEnRWTh7KFqY5Yr-JVZaV-Xje8CXsVRjGHOkt7cBXwQwy7mAYN83rb0j0mJqntgjbs1y4VBVOGbd9XmvoRe0iJkN9H9IUGXMbadAsFKXJKz83DxcDAagCLeaaFJljx-5lysMX2iarwpWDk66k6059bmwdc_Cv72X2T4cld8hmIKYNsDCUn-aQ7v9cRN5UrRs49sZaMmrOxp5lEltakeOxBSxm9_hPcpTqPj3OEndA"
      val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")
      JwtJson.decodeJson(auth0Token, auth0PublicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "fail when the JWT is not valid" in {
      cancel("todo")
      //decodeJWT_RSA("what?", "secret?") mustBe a[Failure[_]]
    }

  }

}
