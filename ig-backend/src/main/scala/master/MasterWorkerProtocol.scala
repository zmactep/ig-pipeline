package master

import akka.actor.ActorRef

/**
 * Created with IntelliJ IDEA.
 * User: Kos
 * Date: 25.09.13
 * Time: 11:35
 * To change this template use File | Settings | File Templates.
 */
object MasterWorkerProtocol {
  // Messages from Workers
  case class WorkerCreated(worker: ActorRef)
  case class WorkerRequestsWork(worker: ActorRef)
  case class WorkIsDone(worker: ActorRef)

  // Messages to Workers
  case class WorkToBeDone(work: Any)
  case object WorkIsReady
  case object NoWorkToBeDone
}
