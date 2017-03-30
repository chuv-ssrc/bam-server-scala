package scala.auth

import pdi.jwt.{Jwt, JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import org.scalatestplus.play._
import auth.Auth._
import auth.RSA._
import java.security._
import models.User
import scala.setup.WithToken



class AuthSpec extends PlaySpec with WithToken {

  /* To generate id_rsa, id_rsa.pub test files: */
  //writeTestKeyPair(512)

  val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")

  "JWT decoding" should {

    "return the claim when the JWT is valid" in {
      val token = Jwt.encode("""{"user":1}""", privateKey, JwtAlgorithm.RS256)
      JwtJson.decodeJson(token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "decode an auth0 token" in {
      JwtJson.decodeJson(auth0Token, auth0PublicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

  }

  "validateToken" should {

    "return the claim as a User instance" in {
      val claim = validateToken(auth0Token, "auth0")
      claim mustBe a [Success[_]]
      claim.get mustBe a [User]
    }

    "fail when the JWT is not valid" in {
      cancel("todo")
      //decodeJWT_RSA("what?", "secret?") mustBe a[Failure[_]]
    }

  }

}
