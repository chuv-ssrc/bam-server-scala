package auth

import javax.inject._

import play.api.libs.concurrent.Execution.Implicits._
import play.api.mvc.Results._
import play.api.mvc._
import play.api.http.HeaderNames._

import scala.concurrent.Future
import scala.util.{Failure, Success, Try}
import models.User
import play.api.Logger
import play.api.db.Database


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
@Singleton
class AuthenticatedAction @Inject()(db: Database) extends ActionBuilder[AuthenticatedRequest] {

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
        val Auth = new Auth(db)
        Auth.validateToken(jwt, db) match {
          case Failure(err) =>
            Logger.debug("Failed token validation: "+ err.getMessage)
            Future.successful(Unauthorized(err.getMessage))
          case Success(user: User) =>
            if (!user.isActive) {
              Future.successful(Forbidden("This user requires validation by admin"))
            } else {
              block(new AuthenticatedRequest(user, request))
            }
        }
    }
  }
}

//@Singleton
//class Admin @Inject()(db: Database) extends AuthenticatedAction(db) {
//
//  override def invokeBlock[A](request: Request[A], block: AuthenticatedRequest[A] => Future[Result]) = {
//    super.invokeBlock[A](request, block)
//  }
//
//}

