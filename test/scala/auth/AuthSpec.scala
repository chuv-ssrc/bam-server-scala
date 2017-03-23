package scala.auth

import pdi.jwt.{Jwt, JwtAlgorithm, JwtJson}

import scala.util.{Failure, Success, Try}
import play.api.libs.json._
import org.scalatestplus.play._
import auth.Auth._
import auth.RSA._

import java.security._


class AuthSpec extends PlaySpec {

  val privateKey: PrivateKey = readPrivateKeyFromPemFile("resources/rsa_keys/test/id_rsa")
  val publicKey: PublicKey = readPublicKeyFromPemFile("resources/rsa_keys/test/id_rsa.pub")

  //writeTestKeyPair(512)

  "`validateJWT`" should {

    "return the claim when the JWT is valid" in {
      val token = Jwt.encode("""{"user":1}""", privateKey, JwtAlgorithm.RS256)
      JwtJson.decodeJson(token, publicKey, Seq(JwtAlgorithm.RS256)) mustBe a[Success[_]]
    }

    "fail when the JWT is not valid" in {
      //decodeJWT_RSA("what?", "secret?") mustBe a[Failure[_]]
    }

  }

}
