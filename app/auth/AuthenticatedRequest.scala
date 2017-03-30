package auth

import play.api.libs.concurrent.Execution.Implicits._
import play.api.mvc.Results._
import play.api.mvc._
import play.api.http.HeaderNames._

import scala.concurrent.Future
import scala.util.{Success, Failure, Try}
import models.User
import auth.Auth.validateToken
import play.api.Logger


// N.B. "Unauthorized" (401) means "unauthenticated". "Forbidden" (403) means "not enought rights".

/**
  * Extends the Request object to make it require a parameter of type User.
  * An `AuthenticatedAction` will only execute an `AuthenticatedRequest` after extracting a User
  * from the initial Request.
  */
class AuthenticatedRequest[A](val user: User, request: Request[A]) extends WrappedRequest[A](request)

/**
  * Checks that the JWT in the Authorization header of a POST request contains a known User.
  * If there is a User but he does not have the right to access the resource, return status is Forbidden (403).
  * If there is no User, return status is Unauthorized (401).
  */
object AuthenticatedAction extends ActionBuilder[AuthenticatedRequest] {

  def jwtFromAuthHeader[A](request: Request[A]): Option[String] = {
    request.headers.get(AUTHORIZATION) map (_.split(" ").last)  // remove "Bearer" prefix
  }

  def jwtFromUrlParam[A](request: Request[A]): Option[String] = {
    request.getQueryString("token")
  }

  def invokeBlock[A](request: Request[A], block: AuthenticatedRequest[A] => Future[Result]) = {
    val maybeJwt: Option[String] = jwtFromAuthHeader(request) orElse jwtFromUrlParam(request)
    maybeJwt match {
      case None =>
        Logger.debug("No Authorization header")
        Future.successful(Unauthorized("No Authorization header"))
      case Some(jwt: String) =>
        validateToken(jwt, "auth0") match {
          case Failure(err) =>
            Logger.debug("Failed token validation: "+ err.getMessage)
            Future.successful(Unauthorized(err.getMessage))
          case Success(user: User) =>
            block(new AuthenticatedRequest(user, request))
        }
    }
  }
}