package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.parse.ParseClassName;
import com.parse.ParseObject;

import java.util.Date;

/**
 * Created by ian on 23/07/16.
 */
@ParseClassName(Constants.RUNNER)
public class RunnerDTO extends ParseObject{

    public EntityDTO getJockey(){
        return (EntityDTO)getParseObject(Constants.JOCKEY);
    }
    public EntityDTO getTrainer(){
        return (EntityDTO)getParseObject(Constants.TRAINER);
    }
    public EntityDTO getHorse(){
        return (EntityDTO)getParseObject(Constants.HORSE);
    }
    public String getSp(){
        return getString(Constants.SP);
    }
    public String getOr(){
        return getString(Constants.OR);
    }
    public String getStall(){
        return getString(Constants.STALL);
    }
    public String getSilkImgLink(){
        return getString(Constants.SILK_IMG_LINK);
    }
    public String getWeight(){
        return getString(Constants.WEIGHT);
    }
    public String getBreeding(){
        return getString(Constants.BREEDING);
    }
    public String getForm(){
        return getString(Constants.FORM);
    }
    public String getAge(){
        return getString(Constants.AGE);
    }
    public Boolean isRunning(){
        return getBoolean(Constants.RUNNING);
    }
    public Date getCreatedAt(){
        return getDate(Constants.CREATEDAT);
    }

    public void setJockey(EntityDTO jockey){
        put(Constants.JOCKEY, jockey);
    }
    public void setTrainer(EntityDTO trainer){put(Constants.TRAINER, trainer);}
    public void setHorse(EntityDTO horse){
        put(Constants.HORSE, horse);
    }
    public void setSp(String sp){put(Constants.SP, sp);}
    public void setOr(String or){
        put(Constants.OR, or);
    }
    public void setStall(String stall){
        put(Constants.STALL, stall);
    }
    public void setSilkImgLink(String silkImgLink){
        put(Constants.SILK_IMG_LINK, silkImgLink);
    }
    public void setWeight(String weight){
        put(Constants.WEIGHT, weight);
    }
    public void setBreeding(String breeding){
        put(Constants.BREEDING, breeding);
    }
    public void setForm(String form){
        put(Constants.FORM, form);
    }
    public void setAge(String age){
        put(Constants.AGE, age);
    }
    public void setRunning(Boolean running){
        put(Constants.RUNNING, running);
    }
    public void setCreatedAt(Date createdAt){
        put(Constants.CREATEDAT, createdAt);
    }

}
