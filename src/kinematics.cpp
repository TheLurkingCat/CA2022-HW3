#include "kinematics.h"

#include <algorithm>

#include "utils.h"
void forwardKinematics(const Posture& posture, Bone* bone) {
  // TODO (FK)
  // Same as HW2, but have some minor change
  // Hint:
  //   1. If you don't use `axis` in this function, you can copy-paste your code
  // Note:
  //   1. bone.axis becomes quaternion instead of vector3f
  while (bone != nullptr) {
    if (bone->parent == nullptr) {
      bone->startPosition = bone->endPosition = posture.translations[0];
      bone->rotation = posture.rotations[0];
    } else {
      bone->rotation = bone->parent->rotation * bone->rotationParentCurrent * posture.rotations[bone->idx];
      bone->startPosition = bone->parent->endPosition + posture.translations[bone->idx];
      bone->endPosition = bone->startPosition + bone->rotation * bone->direction * bone->length;
    }
    forwardKinematics(posture, bone->child);
    bone = bone->sibling;
  }
}

Eigen::VectorXf leastSquareSolver(const Eigen::Matrix3Xf& jacobian, const Eigen::Vector3f& target) {
  // TODO (find x which min(| jacobian * x - target |))
  // Hint:
  //   1. Linear algebra - least squares solution
  //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
  // Note:
  //   1. SVD or other pseudo-inverse method is useful
  //   2. Some of them have some limitation, if you use that method you should check it.
  return jacobian.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(target);
}

void inverseKinematics(const Eigen::Vector3f& target, Bone* start, Bone* end, Posture& posture) {
  constexpr int maxIterations = 10000;
  constexpr float epsilon = 1E-3;
  constexpr float step = 0.1;
  // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is root.
  Bone* root = start - start->idx;
  std::vector<Bone*> boneList;
  // TODO
  // Hint:
  //   1. Traverse from end to start is easier than start to end (since there is only 1 parent)
  // Note:
  //   1. Both start and end should be in the list
  //   2. start bone is guaranteed to be reachable from end bone.
  Bone* current = end;
  while (current != nullptr && current != start->parent) {
    boneList.emplace_back(current);
    current = current->parent;
  }
  // End TODO
  size_t boneNum = boneList.size();
  Eigen::Matrix3Xf jacobian(3, 3 * boneNum);
  jacobian.setZero();

  for (int i = 0; i < maxIterations; ++i) {
    forwardKinematics(posture, root);
    Eigen::Vector3f desiredVector = target - end->endPosition;
    if (desiredVector.norm() < epsilon) break;
    // TODO (compute jacobian)
    // Hint:
    //   1. You should not put rotation in jacobian if it doesn't have that DoF.
    for (size_t j = 0, column = 0; j < boneNum; j++, column += 3) {
      const auto& bone = *boneList[j];
      Eigen::Vector3f armVector = end->endPosition - bone.startPosition;
      Eigen::Matrix3f rotation = bone.rotation.toRotationMatrix();
      if (bone.dofrx) jacobian.col(column + 0) = rotation.col(0).normalized().cross(armVector);
      if (bone.dofry) jacobian.col(column + 1) = rotation.col(1).normalized().cross(armVector);
      if (bone.dofrz) jacobian.col(column + 2) = rotation.col(2).normalized().cross(armVector);
    }
    // End TODO
    Eigen::VectorXf deltaTheta = step * leastSquareSolver(jacobian, desiredVector);
    // TODO (update rotation)
    // Hint:
    //   1. Use .toRotationMatrix().eulerAngles(2, 1, 0).reverse() to get euler angles from quaternion
    //   2. It returns radian, not degree.
    //   3. You can ignore rotation limit of the bone.

    for (size_t j = 0, column = 0; j < boneNum; j++, column += 3) {
      const auto& bone = *boneList[j];
      Eigen::Vector3f angle = posture.rotations[bone.idx].toRotationMatrix().eulerAngles(2, 1, 0).reverse();
      angle += deltaTheta.segment<3>(column);
      posture.rotations[bone.idx] = Eigen::AngleAxisf(angle[2], Eigen::Vector3f::UnitZ()) *
                                    Eigen::AngleAxisf(angle[1], Eigen::Vector3f::UnitY()) *
                                    Eigen::AngleAxisf(angle[0], Eigen::Vector3f::UnitX());
    }
  }
}

// Bonus:
// You can implementation another IK solver, in this homework we use inverse-jacobian solver.
