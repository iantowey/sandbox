package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;

import com.lucidlogic.horsetracker.config.Constants;
import com.parse.ParseClassName;
import com.parse.ParseObject;

import java.util.Date;

public class Runner extends BaseObservable{

    private String id;
    private Entity jockey;
    private Entity trainer;
    private Entity horse;
    private String sp;
    private String or;
    private String silkImgLink;
    private String stall;
    private String weight;
    private String breeding;
    private String form;
    private String age;
    private Boolean running;
    private Date createdAt;


    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public Entity getJockey() {
        return jockey;
    }

    public void setJockey(Entity jockey) {
        this.jockey = jockey;
    }

    public Entity getTrainer() {
        return trainer;
    }

    public void setTrainer(Entity trainer) {
        this.trainer = trainer;
    }

    public Entity getHorse() {
        return horse;
    }

    public void setHorse(Entity horse) {
        this.horse = horse;
    }

    public String getSp() {
        return sp;
    }

    public void setSp(String sp) {
        this.sp = sp;
    }

    public String getWeight() {
        return weight;
    }

    public void setWeight(String weight) {
        this.weight = weight;
    }

    public String getBreeding() {
        return breeding;
    }

    public void setBreeding(String breeding) {
        this.breeding = breeding;
    }

    public String getForm() {
        return form;
    }

    public void setForm(String form) {
        this.form = form;
    }

    public String getAge() {
        return age;
    }

    public void setAge(String age) {
        this.age = age;
    }

    public Boolean getRunning() {
        return running;
    }

    public void setRunning(Boolean running) {
        this.running = running;
    }

    public Date getCreatedAt() {
        return createdAt;
    }

    public void setCreatedAt(Date createdAt) {
        this.createdAt = createdAt;
    }

    public String getOr() {
        return or;
    }

    public void setOr(String or) {
        this.or = or;
    }

    public String getSilkImgLink() {
        return silkImgLink;
    }

    public void setSilkImgLink(String silkImgLink) {
        this.silkImgLink = silkImgLink;
    }

    public String getStall() {
        return stall;
    }

    public void setStall(String stall) {
        this.stall = stall;
    }
}
