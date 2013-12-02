/*

	Copyright 2010 Etay Meiri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
    using namespace std;
#include "pipeline.h"


const Matrix4f* Pipeline::GetTrans()
{
    //cout<<"lol"<<endl;
    Matrix4f ScaleTrans, RotateTrans, TranslationTrans, CameraTranslationTrans, CameraRotateTrans, PersProjTrans;
   // cout<<m_scale.x<<" "<<m_scale.y<<" "<<m_scale.z<<endl;
    ScaleTrans.InitScaleTransform(m_scale.x, m_scale.y, m_scale.z);
   // cout<<m_rotateInfo.x<<" "<<m_rotateInfo.y<<" "<<m_rotateInfo.z<<endl;
    RotateTrans.InitRotateTransform(m_rotateInfo.x, m_rotateInfo.y, m_rotateInfo.z);
   // cout<<m_worldPos.x<<" "<<m_worldPos.y<<" "<<m_worldPos.z<<endl;
    TranslationTrans.InitTranslationTransform(m_worldPos.x, m_worldPos.y, m_worldPos.z);
    //cout<<m_camera.Pos.x<<" "<<m_camera.Pos.y<<" "<<m_camera.Pos.z<<endl;
    CameraTranslationTrans.InitTranslationTransform(-m_camera.Pos.x, -m_camera.Pos.y, -m_camera.Pos.z);
  //  cout<<m_camera.Target.x<<" "<<m_camera.Target.y<<" "<<m_camera.Target.z<<endl;
 //   cout<<m_camera.Up.x<<" "<<m_camera.Up.y<<" "<<m_camera.Up.z<<endl;
    CameraRotateTrans.InitCameraTransform(m_camera.Target, m_camera.Up);
  //  cout<<m_persProj.FOV<<" "<<m_persProj.Width<<" "<<m_persProj.Height<<m_persProj.zNear<<" "<<m_persProj.zFar<<endl;
    PersProjTrans.InitPersProjTransform(m_persProj.FOV, m_persProj.Width, m_persProj.Height, m_persProj.zNear, m_persProj.zFar);

 /*  for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<CameraTranslationTrans.m[i][j]<<" ";
        }
        cout<<endl;
    }*/
    m_transformation =  /*PersProjTrans */ CameraTranslationTrans * TranslationTrans * RotateTrans * ScaleTrans;
/*    cout<<"begin"<<endl;
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            cout<<m_transformation.m[i][j]<<" ";
        }
        cout<<endl;
    }*/
    return &m_transformation;
}

const Matrix4f* Pipeline::getPerspective()
{
    Matrix4f PersProjTrans;
     PersProjTrans.InitPersProjTransform(m_persProj.FOV, m_persProj.Width, m_persProj.Height, m_persProj.zNear, m_persProj.zFar);
     m_perspective = PersProjTrans;
     return &m_perspective;
}


