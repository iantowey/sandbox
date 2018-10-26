package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.parse.ParseClassName;
import com.parse.ParseObject;

import java.util.Date;
import java.util.List;
import java.util.Map;

/**
 * Created by itowey on 12/08/16.
 */
@ParseClassName(Constants.ENTITYDETAILFUTURE)
public class EntityDetailFutureDTO extends ParseObject{

    public  String getSex(){return getString(Constants.SEX);}
    public String getAge(){return getString(Constants.AGE);}
    public Map getSire(){return getMap(Constants.SIRE);}
    public Map getDam(){return getMap(Constants.DAM);}
    public String getName(){return getString(Constants.NAME);}
    public String getOwner(){return getString(Constants.OWNER);}
    public List getFutEnt(){return getList(Constants.FUT_ENT);}
    public Date getCollectionDate(){return getDate(Constants.COLLECTION_DATE);}

    public  void setSex(String sex){put(Constants.SEX, sex);}
    public void setAge(String age){put(Constants.AGE, age);}
    public void setSire(Map sire){put(Constants.SIRE, sire);}
    public void setDam(Map dam){put(Constants.DAM, dam);}
    public void setName(String name){put(Constants.NAME, name);}
    public void setOwner(String owner){put(Constants.OWNER, owner);}
    public void setFutEnt(List list){put(Constants.FUT_ENT, list);}
    public void setCollectionDate(Date date){put(Constants.COLLECTION_DATE, date);}

}
